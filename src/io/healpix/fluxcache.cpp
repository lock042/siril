#include "config.h"
#include "fluxcache.h"
#include <cstdint>

#ifdef HAVE_SQLITE
#include <sqlite3.h>
#endif

#include <ctime>
#include <filesystem>
#include <mutex>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <omp.h>
#include "core/siril_log.h"
#include "core/siril.h"

extern "C" {
#include "core/siril_app_dirs.h"
}

std::map<uint32_t, std::weak_ptr<FluxCache>> FluxCache::instances;
std::mutex FluxCache::mtx;

// public utility function to ensure we have our cache directory
std::filesystem::path FluxCache::get_or_create_cache_dir() {
    const char *uddir = siril_get_user_data_dir();
    std::filesystem::path siril_subdir = std::filesystem::path(uddir) / "spcc-cache";
    try {
        if (!std::filesystem::exists(siril_subdir)) {
            std::filesystem::create_directories(siril_subdir);
        }
        return siril_subdir;
    } catch (const std::filesystem::filesystem_error& e) {
        siril_debug_print(_("Can't create directory %s, falling back to system temp\n"),
                          siril_subdir.string().c_str());
        return std::filesystem::temp_directory_path();
    }
}


#ifndef HAVE_SQLITE
/* -------------------------------------------------------------------------
   STUB IMPLEMENTATION (SQLite Disabled)
   ------------------------------------------------------------------------- */

FluxCache::FluxCache(const std::string& dbPath) : db(nullptr) { }
FluxCache::~FluxCache() { }

bool FluxCache::getCacheEntry(uint32_t hp_id, std::vector<char>& out_data) { return false; }
void FluxCache::setCacheEntry(const std::vector<CacheSegment>& segments) { }
void FluxCache::purge(int days) { }
FluxCacheStats FluxCache::getDatabaseStats() { return {0, 0, 0}; }
std::vector<int> FluxCache::getCachedHealpixelIds() { return {} };

std::shared_ptr<FluxCache> FluxCache::getCache(uint32_t c_lvl, uint32_t h_lvl, uint32_t chunk_id) {
    // Return a dummy object so calling code doesn't crash on null
    return std::shared_ptr<FluxCache>(new FluxCache(""));
}

void FluxCache::closeAllCaches() {
    std::lock_guard<std::mutex> lock(mtx);
    instances.clear();
}

#else

// Note, constructor is private. Use the static factory
FluxCache::FluxCache(const std::string& dbPath) {
    bool corrupted = false;
    std::error_code ec;

    if (std::filesystem::exists(dbPath, ec)) {
        // Quick Magic Header Check: SQLite files must start with "SQLite format 3\000"
        std::ifstream f(dbPath, std::ios::binary);
        char magic[16];
        if (!f.read(magic, 16) || std::memcmp(magic, "SQLite format 3\000", 16) != 0) {
            siril_log_color_message(_("FluxCache: Database corrupted. Resetting...\n"), "salmon");
            corrupted = true;
        }
        f.close();

        // Structural Integrity Check
        if (!corrupted) {
            sqlite3* test_db = nullptr;
            if (sqlite3_open_v2(dbPath.c_str(), &test_db, SQLITE_OPEN_READONLY, nullptr) == SQLITE_OK) {
                // quick_check is fast and catches the most common corruption issues
                int rc = sqlite3_exec(test_db, "PRAGMA quick_check;", nullptr, nullptr, nullptr);
                if (rc != SQLITE_OK) {
                    siril_log_color_message(_("FluxCache: Database corrupted. Resetting...\n"), "salmon");
                    corrupted = true;
                }
                sqlite3_close(test_db);
            } else {
                corrupted = true;
            }
        }
    }

    if (corrupted) {
        // Delete the main DB and the WAL/SHM side-car files
        std::filesystem::remove(dbPath, ec);
        std::filesystem::remove(dbPath + "-wal", ec);
        std::filesystem::remove(dbPath + "-shm", ec);
    }

    // Open or create the database
    if (sqlite3_open(dbPath.c_str(), &db) == SQLITE_OK) {
        const char* sql = "CREATE TABLE IF NOT EXISTS flux_cache ("
                          "hp_id INTEGER PRIMARY KEY, "
                          "timestamp INTEGER, "
                          "data BLOB);";
        sqlite3_exec(db, sql, nullptr, nullptr, nullptr);
        
        // Optimize for blob performance
        sqlite3_exec(db, "PRAGMA synchronous = OFF;", nullptr, nullptr, nullptr);
        sqlite3_exec(db, "PRAGMA journal_mode = WAL;", nullptr, nullptr, nullptr);
    } else {
        siril_log_color_message(_("FluxCache: Failed to open/create database at %s\n"), "red", dbPath.c_str());
    }
}


// Close all database handles that are held open in our singleton
void FluxCache::closeAllCaches() {
    std::lock_guard<std::mutex> lock(mtx);
    siril_log_message(_("FluxCache: Closing all cached database handles...\n"));
    instances.clear();
}

// Close database on descruction
FluxCache::~FluxCache() {
    if (db) {
        // Merge WAL data and truncate the WAL file to 0 bytes
        sqlite3_wal_checkpoint_v2(db, NULL, SQLITE_CHECKPOINT_TRUNCATE, NULL, NULL);
        sqlite3_close_v2(db); // Use v2 for safer cleanup
        db = nullptr;
    }
}

// Get a cache instance for the given healpixel parameters
std::shared_ptr<FluxCache> FluxCache::getCache(uint32_t chunk_level, uint32_t healpix_level, uint32_t chunk_id) {
    std::lock_guard<std::mutex> lock(mtx);

    // Try to find the instance
    auto it = instances.find(chunk_id);
    if (it != instances.end()) {
        if (auto shared_cache = it->second.lock()) {
            return shared_cache; 
        }
        // If lock() fails, it means the object was destroyed,
        // so we fall through to create a new one.
    }

    // Not found or expired: Create new instance
    std::filesystem::path cache_dir = get_or_create_cache_dir();
    std::string db_filename = "siril_cat" + std::to_string(chunk_level) + "_healpix" +
                              std::to_string(healpix_level) + "_xpsamp_" +
                              std::to_string(chunk_id) + ".dat.cache";
    std::string full_path = (cache_dir / db_filename).string();

    siril_log_message(_("FluxCache: Initializing database for chunk %u at %s\n"),
                      chunk_id, full_path.c_str());

    auto cache = std::shared_ptr<FluxCache>(new FluxCache(full_path));
    instances[chunk_id] = cache;

    return cache;
}

// Get a cache entry
bool FluxCache::getCacheEntry(uint32_t hp_id, std::vector<char>& out_data) {
    sqlite3_stmt* stmt;
    bool found = false;

    // Get the data
    const char* select_sql = "SELECT data FROM flux_cache WHERE hp_id = ?;";
    if (sqlite3_prepare_v2(db, select_sql, -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_int64(stmt, 1, (sqlite3_int64)hp_id);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            found = true;
            const void* blob = sqlite3_column_blob(stmt, 0);
            int size = sqlite3_column_bytes(stmt, 0);
            if (size > 0 && blob) {
                out_data.assign((const char*)blob, (const char*)blob + size);
            }
        }
        sqlite3_finalize(stmt);
    }

    // Conditional update (Only writes to disk if the day has changed)
    if (found) {
        const char* update_sql = "UPDATE flux_cache SET timestamp = strftime('%s','now') "
                                 "WHERE hp_id = ? AND date(timestamp, 'unixepoch') != date('now');";

        if (sqlite3_prepare_v2(db, update_sql, -1, &stmt, nullptr) == SQLITE_OK) {
            sqlite3_bind_int64(stmt, 1, (sqlite3_int64)hp_id);
            sqlite3_step(stmt);
            sqlite3_finalize(stmt);
        }
    }

    return found;
}

// Set a cache entry
void FluxCache::setCacheEntry(const std::vector<CacheSegment>& segments) {
    if (!db || segments.empty()) return;

    sqlite3_stmt* stmt;
    const char* sql = "INSERT OR REPLACE INTO flux_cache (hp_id, timestamp, data) VALUES (?, ?, ?);";
    
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) return;

    sqlite3_exec(db, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
    long long now = (long long)std::time(nullptr);

    for (const auto& seg : segments) {
        sqlite3_bind_int(stmt, 1, seg.hp_id);
        sqlite3_bind_int64(stmt, 2, now);
        sqlite3_bind_blob(stmt, 3, seg.data, static_cast<int>(seg.size), SQLITE_STATIC);
        
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }

    sqlite3_exec(db, "END TRANSACTION;", nullptr, nullptr, nullptr);
    sqlite3_finalize(stmt);
}

// Private utility function to gather the paths of our cache shards
std::vector<std::filesystem::path> FluxCache::getCachePaths() {
    std::vector<std::filesystem::path> paths;
    std::filesystem::path cache_dir = get_or_create_cache_dir();

    if (!std::filesystem::exists(cache_dir)) return paths;

    for (const auto& entry : std::filesystem::directory_iterator(cache_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".cache") {
            paths.push_back(entry.path());
        }
    }
    return paths;
}

// Purge all cache entries older than a certain number of days
void FluxCache::purge(int days) {
    // if (days < 0) return;
    const long long expiry_seconds = static_cast<long long>(days) * 86400LL;
    std::filesystem::path cache_dir = get_or_create_cache_dir();

    if (!std::filesystem::exists(cache_dir)) return;

    // Initial Stats
    FluxCacheStats stats = getDatabaseStats();
    char *readable_size = g_format_size(stats.total_size_bytes);
    siril_log_message(_("Starting clean up. Total cache files = %d, Total size in bytes = %s, Total number of rows = %d\n"),
                      stats.file_count, readable_size, stats.total_rows);
    g_free(readable_size);

    std::vector<std::filesystem::path> paths = getCachePaths();

    const std::string delete_sql =
        "DELETE FROM flux_cache WHERE timestamp < (strftime('%s','now') - " +
        std::to_string(expiry_seconds) + ");";

    std::vector<std::filesystem::path> paths_to_vacuum;
    int total_deleted = 0;

    // PHASE 1: Fast Delete (All Threads)
    // We scan everything and delete expired rows.
    siril_log_message(_("Removing rows older than %d days\n"), days);
    #pragma omp parallel for reduction(+:total_deleted)
    for (size_t i = 0; i < paths.size(); ++i) {
        sqlite3* temp_db = nullptr;
        if (sqlite3_open(paths[i].string().c_str(), &temp_db) == SQLITE_OK) {
            sqlite3_busy_timeout(temp_db, 1000);
            if (sqlite3_exec(temp_db, delete_sql.c_str(), nullptr, nullptr, nullptr) == SQLITE_OK) {
                int count = sqlite3_changes(temp_db);
                if (count > 0) {
                    total_deleted += count;
                    // Store this path for the second phase
                    #pragma omp critical
                    paths_to_vacuum.push_back(paths[i]);
                }
            }
            sqlite3_close(temp_db);
        }
    }

    // PHASE 2: Heavy Vacuum (2 Threads)
    // Only touch the files that actually changed.
    if (!paths_to_vacuum.empty()) {
        siril_log_message(_("Reclaiming disk space from %zu modified cache files...\n"), paths_to_vacuum.size());

        #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < paths_to_vacuum.size(); ++i) {
            sqlite3* temp_db = nullptr;
            if (sqlite3_open(paths_to_vacuum[i].string().c_str(), &temp_db) == SQLITE_OK) {
                sqlite3_exec(temp_db, "VACUUM;", nullptr, nullptr, nullptr);
                sqlite3_close(temp_db);
            }
        }
    }

    // Final Stats
    stats = getDatabaseStats();
    readable_size = g_format_size(stats.total_size_bytes);
    siril_log_message(_("Clean up complete. Total entries removed: %d, Total cache files = %d, Total size in bytes = %s, Total number of rows = %d\n"), total_deleted, stats.file_count, readable_size, stats.total_rows);
    g_free(readable_size);
}

// Return database statistics
FluxCacheStats FluxCache::getDatabaseStats() {
    uint32_t f_count = 0;
    uint64_t s_bytes = 0;
    uint64_t r_rows = 0;

    std::filesystem::path cache_dir = get_or_create_cache_dir();

    if (!std::filesystem::exists(cache_dir)) return {0, 0, 0};

    std::vector<std::filesystem::path> paths = getCachePaths();

    #pragma omp parallel for reduction(+:f_count, s_bytes, r_rows)
    for (size_t i = 0; i < paths.size(); ++i) {
        f_count++;
        s_bytes += std::filesystem::file_size(paths[i]);

        sqlite3* temp_db = nullptr;
        if (sqlite3_open_v2(paths[i].string().c_str(), &temp_db, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, nullptr) == SQLITE_OK) {
            sqlite3_stmt* stmt;
            if (sqlite3_prepare_v2(temp_db, "SELECT count(*) FROM flux_cache", -1, &stmt, nullptr) == SQLITE_OK) {
                if (sqlite3_step(stmt) == SQLITE_ROW) {
                    r_rows += sqlite3_column_int64(stmt, 0);
                }
                sqlite3_finalize(stmt);
            }
            sqlite3_close(temp_db);
        }
    }

    return { f_count, s_bytes, r_rows };
}

// Return a list of all cached healpixels
std::vector<int> FluxCache::getCachedHealpixelIds() {
    std::vector<int> all_cached_ids;
    std::vector<std::filesystem::path> paths = getCachePaths();

    #pragma omp parallel
    for (const auto& path : paths) {
        sqlite3* temp_db = nullptr;
        if (sqlite3_open(path.string().c_str(), &temp_db) == SQLITE_OK) {
            const char* sql = "SELECT hp_id FROM flux_cache;";
            sqlite3_stmt* stmt;

            if (sqlite3_prepare_v2(temp_db, sql, -1, &stmt, nullptr) == SQLITE_OK) {
                while (sqlite3_step(stmt) == SQLITE_ROW) {
                    int hp_id = sqlite3_column_int(stmt, 0);
                    #pragma omp critical
                    all_cached_ids.push_back(hp_id);
                }
                sqlite3_finalize(stmt);
            }
            sqlite3_close(temp_db);
        }
    }

    std::sort(all_cached_ids.begin(), all_cached_ids.end());
    return all_cached_ids;
}

#endif // HAVE_SQLITE


// C Interface
extern "C" {

    // To close all database handles on SPCC dialog shutdown
    void flux_cache_close_all() {
        FluxCache::closeAllCaches();
    }

    
    // To purge old records
    void flux_cache_purge(int days) {
        FluxCache::purge(days);
    }

    // To get stats
    void flux_cache_get_stats(uint32_t *file_count, uint64_t *total_bytes, uint64_t *total_rows) {
        if (!file_count || !total_bytes || !total_rows) return;

        FluxCacheStats stats = FluxCache::getDatabaseStats();
        *file_count = stats.file_count;
        *total_bytes = stats.total_size_bytes;
        *total_rows = stats.total_rows;
    }

    // To get a list of healpixels in the entire cache
    int* get_cached_l8_healpixels(int *count) {
        if (!count) return nullptr;

        std::vector<int> ids = FluxCache::getCachedHealpixelIds();

        *count = static_cast<int>(ids.size());
        if (ids.empty()) {
            return nullptr;
        }

        int* c_array = static_cast<int*>(g_malloc(ids.size() * sizeof(int)));

        std::copy(ids.begin(), ids.end(), c_array);

        return c_array;
    }

}

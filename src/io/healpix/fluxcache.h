#ifndef FLUX_CACHE_H
#define FLUX_CACHE_H

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <cstdint>
#include <mutex>
#include <filesystem>

// Forward declare the header so we don't need the whole project header here
struct HealpixCatHeader;

struct FluxCacheStats {
    size_t file_count;         // Number of .cache files (max 48)
    size_t total_size_bytes;   // Sum of all file sizes on disk
    size_t total_rows;         // Approximate/Quick count of hp_id entries
};

#ifdef HAVE_SQLITE
struct sqlite3;
#endif

class FluxCache {
public:
    struct CacheSegment {
        uint32_t hp_id;
        const char* data;
        size_t size;
    };

private:
#ifdef HAVE_SQLITE
    sqlite3* db = nullptr;
#else
    void* db = nullptr;
#endif

    explicit FluxCache(const std::string& dbPath);
    static std::map<uint32_t, std::weak_ptr<FluxCache>> instances;
    static std::mutex mtx;

public:
    ~FluxCache();

    // Static Factory: Returns the correct DB instance based on chunk_healpix (0-47)
    // For data storage and retrieval
    static std::shared_ptr<FluxCache> getCache(uint32_t chunk_level, uint32_t healpix_level, uint32_t chunk_id);

    // Management function that affect all shards are static functions
    // We don't push sharding complexity to the caller
    static void closeAllCaches();
    static void purge(int days);
    static FluxCacheStats getDatabaseStats();
    static std::filesystem::path get_or_create_cache_dir();
    static std::vector<int> getCachedHealpixelIds();

    // Getter and Setter
    void setCacheEntry(const std::vector<CacheSegment>& segments);
    bool getCacheEntry(uint32_t hp_id, std::vector<char>& out_data);

private:
    static std::vector<std::filesystem::path> getCachePaths();
};

#endif

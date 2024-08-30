#if defined(__APPLE__) && defined(__MACH__)
#import <Foundation/Foundation.h>
#import <sys/mount.h>
#include "FileSystemInfo.h"

int64_t macos_find_space(const char *path) {
    @autoreleasepool {
        NSString *nsPath = [NSString stringWithUTF8String:path];
        struct statfs st;
        int64_t freeSpace;

        if (statfs([nsPath UTF8String], &st) != 0) {
            return -1;
        }

        freeSpace = (int64_t)st.f_bavail * st.f_bsize;

        NSError *error = nil;
        NSDictionary *attributes = [[NSFileManager defaultManager] attributesOfFileSystemForPath:nsPath error:&error];

        if (attributes && !error) {
            NSNumber *purgeableSpace = attributes[NSFileSystemPurgeableSize];
            freeSpace += [purgeableSpace longLongValue];
        }

        return freeSpace;
    }
}
#endif // defined(__APPLE__) && defined(__MACH__)

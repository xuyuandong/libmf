// Ported from chromium.

#ifndef FILE_FILE_ENUMERATOR_H_
#define FILE_FILE_ENUMERATOR_H_

#include <fts.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <stack>
#include <string>
#include <vector>

#include "base/basictypes.h"

namespace file {

// A class for enumerating the files in a provided path. The order of the
// results is not guaranteed.
// This class is blocking and slow, use it carefully.
class FileEnumerator {
 public:
  struct FindInfo {
    struct stat stat;
    std::string filename;
  };

  enum FILE_TYPE {
    FILES                 = 1 << 0,
    DIRECTORIES           = 1 << 1,
  };

  // |root_path| is the starting directory to search for. It may or may not end
  // in a slash.
  //
  // If |recursive| is true, this will enumerate all matches in any
  // subdirectories matched as well. It does a breadth-first search, so all
  // files in one directory will be returned before any files in a
  // subdirectory.
  //
  // |file_type| specifies whether the enumerator should match files,
  // directories, or both.
  //
  // |pattern| is an optional pattern for which files to match. This
  // works like shell globbing. For example, "*.txt" or "Foo???.doc".
  // However, be careful in specifying patterns that aren't cross platform
  // since the underlying code uses OS-specific matching routines.  In general,
  // Windows matching is less featureful than others, so test there first.
  // If unspecified, this will match all files.
  // NOTE: the pattern only matches the contents of root_path, not files in
  // recursive subdirectories.
  // TODO(erikkay): Fix the pattern matching to work at all levels.
  FileEnumerator(const std::string &root_path,
                 bool recursive,
                 FileEnumerator::FILE_TYPE file_type);
  FileEnumerator(const std::string &root_path,
                 bool recursive,
                 FileEnumerator::FILE_TYPE file_type,
                 const std::string &pattern);
  ~FileEnumerator();

  // Returns an empty string if there are no more results.
  const std::string Next();

  // Write the file info into |info|.
  void GetFindInfo(FindInfo* info);

 private:
  std::string root_path_;
  bool recursive_;
  FILE_TYPE file_type_;
  std::string pattern_;  // Empty when we want to find everything.

  // Set to true when there is a find operation open. This way, we can lazily
  // start the operations when the caller calls Next().
  bool is_in_find_op_;

  // A stack that keeps track of which subdirectories we still need to
  // enumerate in the breadth-first search.
  std::stack<std::string> pending_paths_;

  // Returns true if the given path should be skipped in enumeration.
  bool ShouldSkip(const std::string& path);

  FTS* fts_;
  FTSENT* fts_ent_;

  DISALLOW_COPY_AND_ASSIGN(FileEnumerator);
};

// utility method
void MatchFile(const std::string &patterns, std::vector<std::string> *files);

}  // namespace file

#endif  // FILE_FILE_ENUMERATOR_H_

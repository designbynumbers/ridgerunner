/*
 * portability.c — small cross-platform helpers for ridgerunner.
 */

#include "portability.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
#  include <windows.h>
#endif

static void die_io(const char *what, const char *path, const char *file, int line)
{
  fprintf(stderr, "%s (%d): %s '%s' failed: %s\n",
          file, line, what, path, strerror(errno));
  exit(1);
}

void rr_mkdir_or_die(const char *path, const char *file, int line)
{
  int rc;
#if defined(_WIN32) || defined(__MINGW32__)
  rc = mkdir(path);
#else
  rc = mkdir(path, 0755);
#endif
  if (rc != 0 && errno != EEXIST) {
    die_io("mkdir", path, file, line);
  }
}

static void rmtree_contents(const char *path, const char *file, int line)
{
  DIR *d = opendir(path);
  struct dirent *ent;
  char child[4096];
  struct stat st;

  if (d == NULL) {
    return;
  }

  while ((ent = readdir(d)) != NULL) {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
      continue;
    }
    snprintf(child, sizeof(child), "%s/%s", path, ent->d_name);
    if (stat(child, &st) != 0) {
      continue;
    }
    if (S_ISDIR(st.st_mode)) {
      rmtree_contents(child, file, line);
      if (rmdir(child) != 0 && errno != ENOENT) {
        closedir(d);
        die_io("rmdir", child, file, line);
      }
    } else {
      if (remove(child) != 0 && errno != ENOENT) {
        closedir(d);
        die_io("remove", child, file, line);
      }
    }
  }
  closedir(d);
}

void rr_rmtree_or_die(const char *path, const char *file, int line)
{
  struct stat st;
  if (stat(path, &st) != 0) {
    if (errno == ENOENT) {
      return;
    }
    die_io("stat", path, file, line);
  }
  if (S_ISDIR(st.st_mode)) {
    rmtree_contents(path, file, line);
    if (rmdir(path) != 0 && errno != ENOENT) {
      die_io("rmdir", path, file, line);
    }
  } else {
    if (remove(path) != 0 && errno != ENOENT) {
      die_io("remove", path, file, line);
    }
  }
}

int rr_setenv(const char *name, const char *value, int overwrite)
{
#if defined(HAVE_SETENV) && !defined(_WIN32)
  return setenv(name, value, overwrite);
#else
  char *buf;
  size_t n;
  if (!overwrite && getenv(name) != NULL) {
    return 0;
  }
  n = strlen(name) + strlen(value) + 2;
  buf = malloc(n);
  if (buf == NULL) {
    return -1;
  }
  snprintf(buf, n, "%s=%s", name, value);
#  if defined(_WIN32)
  return _putenv(buf);
#  else
  return putenv(buf);
#  endif
#endif
}

void rr_temp_path(char *buf, size_t buflen, const char *filename)
{
  const char *tmp = NULL;
#ifdef _WIN32
  char win_tmp[MAX_PATH];
  DWORD n = GetTempPathA(sizeof(win_tmp), win_tmp);
  if (n > 0 && n < sizeof(win_tmp)) {
    tmp = win_tmp;
  }
#endif
  if (tmp == NULL) {
    tmp = getenv("TMPDIR");
  }
  if (tmp == NULL) {
    tmp = getenv("TEMP");
  }
  if (tmp == NULL) {
    tmp = getenv("TMP");
  }
  if (tmp == NULL) {
    tmp = "/tmp";
  }
  snprintf(buf, buflen, "%s%s%s", tmp,
           (tmp[strlen(tmp) - 1] == '/' || tmp[strlen(tmp) - 1] == '\\') ? "" : "/",
           filename);
}

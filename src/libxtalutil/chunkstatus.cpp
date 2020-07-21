#include "chunkstatus.h"
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

static const char *getTmpDir (void) {
  char *tmpdir = getenv("TMPDIR");
  return tmpdir != NULL ? tmpdir  : "/tmp";
}

/*----------------------------------------------------------------------------*/

ChunkStatus::ChunkStatus() : m_status(NOT_INITIALIZED) {
}

/*----------------------------------------------------------------------------*/

ChunkStatus::Status ChunkStatus::init(const char* id, const char* filename, uint32_t num_chunks) {
  m_id = id;
  m_fn = filename;
  m_compl_inf.resize(num_chunks, 0);
  m_status = filename == NULL ? UNUSED : NOT_LOADED;
  m_num_completed = 0;
  return m_status;
}

/*----------------------------------------------------------------------------*/

ChunkStatus::Status ChunkStatus::load() {
  if (m_status == UNUSED) {
    return m_status;
  }
  if (m_status != NOT_LOADED) {
    m_status = ERROR_LOAD;
    return m_status;
  }
  FILE* f = fopen(m_fn.c_str(), "r");
  /* file not existing (yet) */
  if (!f) {
    m_status = OK;
    return m_status;
  }
  char buffer[1024];
  const uint32_t num_header_lines = 1;
  uint32_t line_count = 0;
  /* compare ids */
  if (fgets(buffer, sizeof(buffer), f) != NULL) {
    /* remove newline */
    buffer[strcspn(buffer, "\n")] = 0;
    if (m_id != buffer) {
      fclose(f);
      m_status = ERROR_LOAD_ID;
      return m_status;
    }
    line_count++;
  }
  while (fgets(buffer, sizeof(buffer), f) != NULL) {
    /* check linecount */
    if (line_count - num_header_lines >= m_compl_inf.size()) {
      fclose(f);
      m_status = ERROR_LOAD_NUM_CHUNKS;
      return m_status;
    }
    /* remove newline */
    buffer[strcspn(buffer, "\n")] = 0;
    /* parse line */
    char s[5];
    int ret = sscanf(buffer, "%*u %4s", s);
    //printf("%u -- %d %s\n", line_count, ret, s);
    if (ret == 1) {
      const uint8_t done = strcmp(s, "done") == 0 ? 1 : 0;
      m_compl_inf[line_count-num_header_lines] = done;
      m_num_completed += done;
    }
    line_count++;
  }
  m_status = (m_id != buffer) ? ERROR_LOAD_ID : m_status;
  m_status = OK;
  fclose(f);
  return m_status;
}

/*----------------------------------------------------------------------------*/

bool ChunkStatus::isCompleted(uint32_t id) const {
  if (m_status != OK) {
    return false;
  }
  return m_compl_inf[id];
}

/*----------------------------------------------------------------------------*/

ChunkStatus::Status ChunkStatus::setCompleted(uint32_t id) {
  while (true) { // single loop
    if (m_status != OK) {
      break;
    }
    m_compl_inf[id] = 1;
    m_num_completed++;
    char *tmpfn = NULL;
    {
      std::string t = getTmpDir();
      t += "/xtal_chunkstatus_XXXXXX";
      tmpfn = strdup(t.c_str());
    }
    int fd = mkstemp(tmpfn);
    FILE *f = fdopen(fd, "wx");
    if (!f) {
      break;
    }
    /* write header */
    fprintf(f, "%s\n", m_id.c_str());
    /* write status */
    for (size_t i = 0; i < m_compl_inf.size(); i++) {
      fprintf(f, "%u %s\n", (uint32_t)i, m_compl_inf[i] ? "done" : "-");
    }
    fclose(f);
    rename(tmpfn, m_fn.c_str());
    free(tmpfn);
    break;
  }
  return m_status;
}

float ChunkStatus::getProgress() const {
  return(float(m_num_completed) / m_compl_inf.size());
}


#ifdef HAS_MAIN

// g++ -D HAS_MAIN chunkstatus.cpp -o /tmp/cs && valgrind --tool=memcheck /tmp/cs && cat -A /tmp/cs.txt

int main(int argc, char const *argv[]) {
  ChunkStatus cs;
  
  const uint32_t num_chunks = 5;
  ChunkStatus::Status r;
  r = cs.init("seqid", "/tmp/cs.txt", num_chunks);
  printf("init: %d\n", r);
  r = cs.load();
  printf("status: ");
  for (uint32_t i = 0; i < num_chunks; i++) {
    printf(" %s", cs.isCompleted(i) ? "X" : ".");
  }
  printf("\n");
  r = cs.setCompleted(3);
  printf("set status: %d\n", r);
  r = cs.setCompleted(0);
  printf("set status: %d\n", r);
  r = cs.setCompleted(num_chunks - 1);
  printf("set status: %d\n", r);

  return 0;
}
#endif

#include <stdint.h>
#include <string>
#include <vector>

/*----------------------------------------------------------------------------*/

class ChunkStatus {
public:
  enum Status {
    OK,
    UNUSED,
    NOT_INITIALIZED,
    NOT_LOADED,
    ERROR_LOAD,
    ERROR_LOAD_FILE,
    ERROR_LOAD_ID,
    ERROR_LOAD_NUM_CHUNKS
  };

  ChunkStatus();
  Status init(const char* id, const char* filename, uint32_t num_chunks);
  Status load();
  bool isCompleted(uint32_t id) const;
  Status setCompleted(uint32_t id); // to be called thread-safe
  float getProgress() const;
private:
  std::string m_id;
  std::string m_fn;
  std::vector<char> m_compl_inf;
  Status m_status;
  uint32_t m_num_completed;
};

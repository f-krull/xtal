
#include <stdio.h>
#include <stdlib.h>

#ifdef EXECNAME
int main(int argc, char **argv) {
   EXECNAME *x = new EXECNAME();
   if (start(x, argc, argv) != 0) {
      fprintf(stderr, "error\n");
      delete x;
      return EXIT_FAILURE;
   }
   delete x;
   fprintf(stderr, "exit\n");
   return EXIT_SUCCESS;
}
#endif


#ifndef EXECTEMPLATE_H_
#define EXECTEMPLATE_H_

class ExecTemplate {
public:
   virtual ~ExecTemplate() {};
   virtual int start() = 0;
   virtual const char* getName() const = 0;
   virtual void registerStuff();
   virtual int start(int argc, char **argv);
};

int start(ExecTemplate *e, int argc, char **argv);


#endif /* EXECTEMPLATE_H_ */

#ifndef ILOGOBSERVER_H_
#define ILOGOBSERVER_H_


class ILogObserver {
public:
   virtual void getLog(const char* msg) = 0;
};

#endif /* ILOGOBSERVER_H_ */

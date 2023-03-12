#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstdarg>

class Gnuplot_Base
{
protected:
  
  std::string exec;
  std::string cmdline;
  
  FILE * fp;
  
public:
  explicit Gnuplot_Base(const char * execname=NULL) : fp(NULL)
  {
    if(execname)
      {
        exec = execname;
      }
    else
      {
          exec="/usr/local/bin/gnuplot";
          //exec = "/sw/bin/gnuplot";
      }
  }
  
  virtual ~Gnuplot_Base()
  {
    end();
  }

  enum { GNUPLOT_OK = 0,
         GNUPLOT_NG,
         GNUPLOT_ALREADY_OPEN,
         GNUPLOT_NOT_OPEN,
         GNUPLOT_CANNOT_EXEC
  };
  
  virtual int begin(const char * cmd_opt=NULL, const char * reserved=NULL) = 0;
  virtual int end() {return 0;}

  int flush()
  {
    if( fp ) ::fflush(fp);
    
    return GNUPLOT_OK;
  }

  int begin_data()
  {
    if( ! fp ) return GNUPLOT_NOT_OPEN;
    return flush();
  }

  int end_data()
  {
    if( ! fp ) return GNUPLOT_NOT_OPEN;
    ::fprintf(fp, "e\n");
    return flush();
  }
  
  int command(const char * cmd, ...) 
  {
    if(!fp) return GNUPLOT_NOT_OPEN;
    if(!cmd) return GNUPLOT_OK;

    { va_list ap; va_start(ap, cmd); ::vfprintf(fp, cmd, ap); va_end(ap); }
    
    return GNUPLOT_OK;
  }

  int commandln(const char * cmd, ...)
  {
    if(!fp) return GNUPLOT_NOT_OPEN;
    if(!cmd) return GNUPLOT_OK;
    
    { va_list ap; va_start(ap, cmd); ::vfprintf(fp, cmd, ap); va_end(ap); }
    ::fprintf(fp, "\n");
    
    return GNUPLOT_OK;
  }
};

class Gnuplot_Pipe : public Gnuplot_Base
{
public:
  virtual int begin(const char * cmd_opt=NULL, const char * reserved=NULL)
  {
    if(fp) return GNUPLOT_ALREADY_OPEN;

    cmdline = exec;
    if(cmd_opt)
      {
        cmdline += " ";
        cmdline += cmd_opt;
      }
    
    fp = ::popen( cmdline.c_str(), "w" );
    if( !fp ) return GNUPLOT_CANNOT_EXEC;

    return GNUPLOT_OK;
  }
  
  virtual int end()
  {
    if( ! fp ) return GNUPLOT_OK;
    
    ::fflush(fp);
    ::pclose(fp);
    fp = NULL;
    
    return GNUPLOT_OK;
  }
  
};


class Gnuplot_Tmpfile : public Gnuplot_Base
{
protected:
  std::string tmpfilename;
  bool remove_file;
  
public:
  virtual int begin(const char * cmd_opt=NULL, const char * filename=NULL)
  {
    if(fp) return GNUPLOT_ALREADY_OPEN;

    cmdline = exec;
    if(cmd_opt)
      {
        cmdline += " ";
        cmdline += cmd_opt;
      }

    if(filename)
      {
        fp = ::fopen(filename, "w");
        tmpfilename = filename;
        remove_file = false;
      }
    else
      {
        char buf[] = "/tmp/gnuplot_XXXXXX";
        int fd = ::mkstemp(buf);
        if( fd == -1 ) return GNUPLOT_NG;
        fp = ::fdopen( fd, "w" );
        tmpfilename = buf;
        remove_file = true;
      }
    
    if( !fp )
      {
        tmpfilename.clear();
        return GNUPLOT_CANNOT_EXEC;
      }
    
    return GNUPLOT_OK;
  }
  
  virtual int end()
  {
    if( ! fp ) return GNUPLOT_OK;
    
    ::fflush(fp);
    ::fclose(fp);
    fp = NULL;
    
    ::system( (cmdline + " " + tmpfilename).c_str() );
    if( remove_file ) ::remove( tmpfilename.c_str() );
    
    return GNUPLOT_OK;
  }
};

// use pipe version as default
typedef Gnuplot_Pipe Gnuplot;

#endif //GNUPLOT_H

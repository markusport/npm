/*! \file cmd_line.h
* \brief simple command line parser
* \author Hanno Hildenbrandt
*/

#ifndef CMD_LINE_H_INCLUDED
#define CMD_LINE_H_INCLUDED

#include <stdexcept>
#include <cassert>
#include <string.h>
#include <sstream>
#include <utility>
#include <string>
#include <type_traits>
#include <regex>
#include "npm.h"


//! \brief namespace of the command line parser
namespace cmd {


  //! \brief Exception type of the command line parser
  class parse_error : public std::invalid_argument
  {
  public:
    parse_error(const char* msg) : std::invalid_argument(msg)
    {}
  };


  //! \brief command line parser
  class cmd_line_parser
  {
  public:
    //! ctor
    cmd_line_parser(int argc, const char** argv)
      : argc_(argc), argv_(argv)
    {}


    //! Returns true if \p name exists in argument, false otherwise
    //! \param name the name of the flag
    bool flag(const char* name) const;


    //! \brief Parse name value pair
    //! \tparam T the type to parse
    //! \param name the name
    //! \param val the value
    //! \param delim delimiter
    //! \returns true if the name-value could be read.
    //!
    //! \p val contains the parsed value on success, otherwise
    //! its unchanged.
    template <typename T>
    bool optional(const char* name, T& val, char delim = '=') const;


    //! \brief Parse name-value pair
    //! \tparam T the return type
    //! \param name the name
    //! \param delim delimiter
    //! \returns the parsed value on success.
    //!
    //! Throws parse_error on failure.
    template <typename T>
    T required(const char* name, char delim = '=') const;

  private:
    int argc_;
    const char** argv_;
  };


  // split argument at delim
  inline std::pair<std::string, std::string> split_arg(const char* carg, char delim)
  {
    const char* s = strchr(carg, delim);
    if (nullptr == s)
    {
      return{ "", "" };
    }
    return{{carg, s}, {s + 1}};
  }


  template <typename T>
  inline void convert_arg(std::pair<std::string, std::string> const& arg, T& x)
  {
    std::istringstream iss(arg.second);
    if (!(iss >> x))
    {
      throw parse_error((std::string("invalid value for argument ") + arg.first).c_str());
    }
  }


  template <>
  inline void convert_arg<npm::Alleles>(std::pair<std::string, std::string> const& arg, npm::Alleles& x)
  {
    std::regex ralleles(R"(\s*(\s*[[:digit:]]+(\.(([[:digit:]]+)?))?)+\s*,?)");
    if (! std::regex_match(arg.second, ralleles))
    {
      throw parse_error((std::string("invalid value for argument ") + arg.first).c_str());
    }
    std::istringstream iss(arg.second);
    for (size_t i = 0; i < npm::Loci::MAX_ALLELE; ++i) 
    {
      iss >> x[i];
      if (!iss)
      {
        throw parse_error((std::string("invalid value for argument ") + arg.first).c_str());
      }
    }
  }


  inline void parse_cmd_flag(const char* name, bool& val, int argc, const char** argv)
  {
    for (int i = 0; i < argc; ++i)
    {
      if (0 == strcmp(argv[i], name)) { val = true; return; }
    }
  }


  template <typename T>
  inline bool parse_optional_arg(const char* name, T& val, int argc, const char** argv, char delim = '=')
  {
    for (int i = 0; i < argc; ++i)
    {
      auto arg = split_arg(argv[i], delim);
      if (arg.first == name) { convert_arg(arg, val); return true; }
    }
    return false;
  }


  template <typename T>
  inline void parse_required_arg(const char* name, T& val, int argc, const char** argv, char delim = '=')
  {
    int i = 0;
    for (; i < argc; ++i)
    {
      auto arg = split_arg(argv[i], delim);
      if (arg.first == name) { convert_arg(arg, val); return; }
    }
    throw parse_error(((std::string("missing argument '") + name) + '\'').c_str());
  }


  inline bool cmd_line_parser::flag(const char* name) const
  {
    bool flag = false;
    parse_cmd_flag(name, flag, argc_, argv_);
    return flag;
  }


  template <typename T>
  inline bool cmd_line_parser::optional(const char* name, T& val, char delim) const
  {
    return parse_optional_arg(name, val, argc_, argv_, delim);
  }


  template <typename T>
  inline T cmd_line_parser::required(const char* name, char delim) const
  {
    T val;
    parse_required_arg(name, val, argc_, argv_, delim);
    return val;
  }


  template <typename N>
  inline int check_any(std::string& p, N const& name, const char* errmsg)
  {
    for (size_t i = 0; i < std::extent<N>::value; ++i)
    {
      if (p == name[i]) return static_cast<int>(i);
    }
    throw parse_error(errmsg);
    return -1;
  }

}

#endif

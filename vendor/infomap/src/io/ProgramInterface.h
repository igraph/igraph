/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef PROGRAM_INTERFACE_H_
#define PROGRAM_INTERFACE_H_

#include "../utils/convert.h"

#include <stdexcept>
#include <utility>
#include <vector>
#include <deque>
#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <unordered_map>

namespace infomap {

class InterruptException : public std::exception {
public:
    InterruptException() {};
};

typedef bool interruptionHandlerFn(void);

struct ArgType {
  static const std::string integer;
  static const std::string number;
  static const std::string string;
  static const std::string path;
  static const std::string probability;
  static const std::string option;
  static const std::string list;
  static const std::unordered_map<std::string, char> toShort;
};

struct Option {
  Option(char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced, bool requireArgument = false, std::string argName = "")
      : shortName(shortName),
        longName(std::move(longName)),
        description(std::move(desc)),
        group(std::move(group)),
        isAdvanced(isAdvanced),
        requireArgument(requireArgument),
        incrementalArgument(false),
        argumentName(std::move(argName)) { }

  virtual ~Option() = default;

  Option(const Option&) = default;
  Option& operator=(const Option&) = default;
  Option(Option&&) = default;
  Option& operator=(Option&&) = default;

  virtual bool parse(std::string const&)
  {
    used = true;
    return true;
  }

  virtual void set(bool value)
  {
    used = true;
    negated = !value;
  };

  Option& setHidden(bool value)
  {
    hidden = value;
    return *this;
  }

  virtual std::ostream& printValue(std::ostream& out) const { return out; }
  virtual std::string printValue() const { return ""; }
  virtual std::string printNumericValue() const { return ""; }

  friend std::ostream& operator<<(std::ostream& out, const Option& option)
  {
    out << option.longName;
    if (option.requireArgument) {
      out << " = ";
      option.printValue(out);
    }
    return out;
  }

  char shortName;
  std::string longName;
  std::string description;
  std::string group;
  bool isAdvanced;
  bool requireArgument;
  bool incrementalArgument;
  std::string argumentName;
  bool hidden = false;
  bool used = false;
  bool negated = false;
};

struct IncrementalOption : Option {
  IncrementalOption(unsigned int& target, char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced)
      : Option(shortName, std::move(longName), std::move(desc), std::move(group), isAdvanced, false), target(target)
  {
    incrementalArgument = true;
  }

  bool parse(std::string const& value) override
  {
    Option::parse(value);
    return ++target;
  }

  void set(bool value) override
  {
    Option::set(value);
    if (value) {
      ++target;
    } else if (target > 0) {
      --target;
    }
  }

  std::ostream& printValue(std::ostream& out) const override { return out << target; }
  std::string printValue() const override { return io::Str() << target; }

  unsigned int& target;
};

template <typename T>
struct ArgumentOption : Option {
  ArgumentOption(T& target, char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced, std::string argName)
      : Option(shortName, std::move(longName), std::move(desc), std::move(group), isAdvanced, true, std::move(argName)), target(target) { }

  bool parse(std::string const& value) override
  {
    Option::parse(value);
    return io::stringToValue(value, target);
  }

  std::string printValue() const override { return io::Str() << target; }
  std::ostream& printValue(std::ostream& out) const override { return out << target; }
  std::string printNumericValue() const override { return TypeInfo<T>::isNumeric() ? printValue() : ""; }

  T& target;
};

template <typename T>
struct LowerBoundArgumentOption : ArgumentOption<T> {
  LowerBoundArgumentOption(T& target, char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced, std::string argName, T minValue)
      : ArgumentOption<T>(target, shortName, std::move(longName), std::move(desc), std::move(group), isAdvanced, std::move(argName)), minValue(minValue) { }

  bool parse(std::string const& value) override
  {
    auto ok = ArgumentOption<T>::parse(value);
    if (ArgumentOption<T>::target < minValue) return false;
    return ok;
  }

  T minValue;
};

template <typename T>
struct LowerUpperBoundArgumentOption : LowerBoundArgumentOption<T> {
  LowerUpperBoundArgumentOption(T& target, char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced, std::string argName, T minValue, T maxValue)
      : LowerBoundArgumentOption<T>(target, shortName, std::move(longName), std::move(desc), std::move(group), isAdvanced, std::move(argName), minValue), maxValue(maxValue) { }

  bool parse(std::string const& value) override
  {
    auto ok = LowerBoundArgumentOption<T>::parse(value);
    if (LowerBoundArgumentOption<T>::target > maxValue) return false;
    return ok;
  }

  T maxValue;
};

template <>
struct ArgumentOption<bool> : Option {
  ArgumentOption(bool& target, char shortName, std::string longName, std::string desc, std::string group, bool isAdvanced)
      : Option(shortName, std::move(longName), std::move(desc), std::move(group), isAdvanced, false), target(target) { }

  bool parse(std::string const& value) override
  {
    Option::parse(value);
    return target = true;
  }

  void set(bool value) override
  {
    Option::set(value);
    target = value;
  }

  std::ostream& printValue(std::ostream& out) const override { return out << target; }
  std::string printValue() const override { return io::Str() << target; }
  std::string printNumericValue() const override { return ""; }

  bool& target;
};

struct ParsedOption {
  explicit ParsedOption(const Option& opt)
      : shortName(opt.shortName),
        longName(opt.longName),
        description(opt.description),
        group(opt.group),
        isAdvanced(opt.isAdvanced),
        requireArgument(opt.requireArgument),
        incrementalArgument(opt.incrementalArgument),
        argumentName(opt.argumentName),
        negated(opt.negated),
        value(opt.printValue()) { }

  friend std::ostream& operator<<(std::ostream& out, const ParsedOption& option)
  {
    if (option.negated)
      out << "no ";
    out << option.longName;
    if (option.requireArgument)
      out << " = " << option.value;
    return out;
  }

  char shortName;
  std::string longName;
  std::string description;
  std::string group;
  bool isAdvanced;
  bool requireArgument;
  bool incrementalArgument;
  std::string argumentName;
  bool negated;
  std::string value;
};

struct TargetBase {
  TargetBase(std::string variableName, std::string desc, std::string group, bool isAdvanced)
      : variableName(std::move(variableName)), description(std::move(desc)), group(std::move(group)), isOptionalVector(false), isAdvanced(isAdvanced) { }

  virtual ~TargetBase() = default;

  TargetBase(const TargetBase&) = default;
  TargetBase& operator=(const TargetBase&) = default;
  TargetBase(TargetBase&&) = default;
  TargetBase& operator=(TargetBase&&) = default;

  virtual bool parse(std::string const& value) = 0;

  std::string variableName;
  std::string description;
  std::string group;
  bool isOptionalVector = false;
  bool isAdvanced;
};

template <typename T>
struct Target : TargetBase {
  Target(T& target, std::string variableName, std::string desc, std::string group, bool isAdvanced)
      : TargetBase(std::move(variableName), std::move(desc), std::move(group), isAdvanced), target(target) { }

  bool parse(std::string const& value) override
  {
    return io::stringToValue(value, target);
  }

  T& target;
};

template <typename T>
struct OptionalTargets : TargetBase {
  OptionalTargets(std::vector<T>& target, std::string variableName, std::string desc, std::string group, bool isAdvanced)
      : TargetBase(std::move(variableName), std::move(desc), std::move(group), isAdvanced), targets(target)
  {
    isOptionalVector = true;
  }

  bool parse(std::string const& value) override
  {
    T target;
    bool ok = io::stringToValue(value, target);
    if (ok)
      targets.push_back(target);
    return ok;
  }

  std::vector<T>& targets;
};

class ProgramInterface {
public:
  ProgramInterface(std::string name, std::string shortDescription, std::string version);

  void setGroups(std::vector<std::string> groups) { m_groups = std::move(groups); }

  template <typename T>
  void addNonOptionArgument(T& target, std::string variableName, std::string desc, std::string group, bool isAdvanced = false)
  {
    m_nonOptionArguments.emplace_back(new Target<T>(target, std::move(variableName), std::move(desc), std::move(group), isAdvanced));
  }

  template <typename T>
  void addOptionalNonOptionArguments(std::vector<T>& target, std::string variableName, std::string desc, std::string group, bool isAdvanced = false)
  {
    if (m_numOptionalNonOptionArguments != 0)
      throw std::runtime_error("Can't have two non-option vector arguments");
    ++m_numOptionalNonOptionArguments;
    m_nonOptionArguments.emplace_back(new OptionalTargets<T>(target, std::move(variableName), std::move(desc), std::move(group), isAdvanced));
  }

  Option& addOptionArgument(char shortName, std::string longName, std::string description, std::string group, bool isAdvanced = false)
  {
    auto* o = new Option(shortName, std::move(longName), std::move(description), std::move(group), isAdvanced);
    m_optionArguments.emplace_back(o);
    return *o;
  }

  Option& addIncrementalOptionArgument(unsigned int& target, char shortName, std::string longName, std::string description, std::string group, bool isAdvanced = false)
  {
    auto* o = new IncrementalOption(target, shortName, std::move(longName), std::move(description), std::move(group), isAdvanced);
    m_optionArguments.emplace_back(o);
    return *o;
  }

  Option& addOptionArgument(bool& target, char shortName, std::string longName, std::string description, std::string group, bool isAdvanced = false)
  {
    auto* o = new ArgumentOption<bool>(target, shortName, std::move(longName), std::move(description), std::move(group), isAdvanced);
    m_optionArguments.emplace_back(o);
    return *o;
  }

  // Without shortName
  Option& addOptionArgument(bool& target, std::string longName, std::string description, std::string group, bool isAdvanced = false)
  {
    return addOptionArgument(target, '\0', std::move(longName), std::move(description), std::move(group), isAdvanced);
  }

  template <typename T>
  Option& addOptionArgument(T& target, char shortName, std::string longName, std::string description, std::string argumentName, std::string group, bool isAdvanced = false)
  {
    auto* o = new ArgumentOption<T>(target, shortName, std::move(longName), std::move(description), std::move(group), isAdvanced, std::move(argumentName));
    m_optionArguments.emplace_back(o);
    return *o;
  }

  template <typename T>
  Option& addOptionArgument(T& target, char shortName, std::string longName, std::string description, std::string argumentName, std::string group, T minValue, bool isAdvanced = false)
  {
    auto* o = new LowerBoundArgumentOption<T>(target, shortName, std::move(longName), std::move(description), std::move(group), isAdvanced, std::move(argumentName), minValue);
    m_optionArguments.emplace_back(o);
    return *o;
  }

  template <typename T>
  Option& addOptionArgument(T& target, char shortName, std::string longName, std::string description, std::string argumentName, std::string group, T minValue, T maxValue, bool isAdvanced = false)
  {
    auto* o = new LowerUpperBoundArgumentOption<T>(target, shortName, std::move(longName), std::move(description), std::move(group), isAdvanced, std::move(argumentName), minValue, maxValue);
    m_optionArguments.emplace_back(o);
    return *o;
  }

  // Without shortName
  template <typename T>
  Option& addOptionArgument(T& target, std::string longName, std::string description, std::string argumentName, std::string group, bool isAdvanced = false)
  {
    return addOptionArgument(target, '\0', std::move(longName), std::move(description), std::move(argumentName), std::move(group), isAdvanced);
  }

  template <typename T>
  Option& addOptionArgument(T& target, std::string longName, std::string description, std::string argumentName, std::string group, T minValue, bool isAdvanced = false)
  {
    return addOptionArgument(target, '\0', std::move(longName), std::move(description), std::move(argumentName), std::move(group), minValue, isAdvanced);
  }

  template <typename T>
  Option& addOptionArgument(T& target, std::string longName, std::string description, std::string argumentName, std::string group, T minValue, T maxValue, bool isAdvanced = false)
  {
    return addOptionArgument(target, '\0', std::move(longName), std::move(description), std::move(argumentName), std::move(group), minValue, maxValue, isAdvanced);
  }

  void parseArgs(const std::string& args);

  std::vector<ParsedOption> getUsedOptionArguments() const;

  unsigned int numRequiredArguments() const { return m_nonOptionArguments.size() - m_numOptionalNonOptionArguments; }

private:
  void exitWithUsage(bool showAdvanced) const;
  void exitWithVersionInformation() const;
  void exitWithError(const std::string& message) const;
  void exitWithJsonParameters() const;

  std::deque<std::unique_ptr<Option>> m_optionArguments;
  std::deque<std::unique_ptr<TargetBase>> m_nonOptionArguments;
  std::string m_programName = "Infomap";
  std::string m_shortProgramDescription;
  std::string m_programVersion;
  std::string m_programDescription;
  std::vector<std::string> m_groups;
  std::string m_executableName = "Infomap";
  unsigned int m_displayHelp = 0;
  bool m_displayVersion = false;
  bool m_printJsonParameters = false;

  unsigned int m_numOptionalNonOptionArguments = 0;
};

} // namespace infomap

#endif // PROGRAM_INTERFACE_H_

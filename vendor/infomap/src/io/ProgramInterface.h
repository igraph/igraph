/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef PROGRAM_INTERFACE_H_
#define PROGRAM_INTERFACE_H_

#include "ParsedOption.h"
#include "../utils/convert.h"

#include <stdexcept>
#include <utility>
#include <vector>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <ostream>
#include <unordered_map>

namespace infomap {

/// Thrown by ProgramInterface help/version/json-parameters paths
/// instead of std::exit(0). The CLI binary catches this and exits with
/// status 0; library callers (R/Python/JS bindings) never trigger these
/// flags. Not derived from std::exception so the generic catch in
/// infomap::run() does not treat it as an error.
struct CleanExit {};

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
  Option(char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_, bool requireArgument_ = false, std::string argName = "")
      : shortName(shortName_),
        longName(std::move(longName_)),
        description(std::move(desc)),
        group(std::move(group_)),
        isAdvanced(isAdvanced_),
        requireArgument(requireArgument_),
        incrementalArgument(false),
        argumentName(std::move(argName)) {}

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

  Option& setChoices(std::vector<std::string> values)
  {
    choices = std::move(values);
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
  std::vector<std::string> choices;
};

struct IncrementalOption : Option {
  IncrementalOption(unsigned int& target_, char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_)
      : Option(shortName_, std::move(longName_), std::move(desc), std::move(group_), isAdvanced_, false), target(target_)
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
  std::string printValue() const override { return io::stringify(target); }

  unsigned int& target;
};

template <typename T>
struct ArgumentOption : Option {
  ArgumentOption(T& target_, char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_, std::string argName)
      : Option(shortName_, std::move(longName_), std::move(desc), std::move(group_), isAdvanced_, true, std::move(argName)), target(target_) {}

  bool parse(std::string const& value) override
  {
    Option::parse(value);
    return io::stringToValue(value, target);
  }

  std::string printValue() const override { return io::stringify(target); }
  std::ostream& printValue(std::ostream& out) const override { return out << target; }
  std::string printNumericValue() const override { return TypeInfo<T>::isNumeric() ? printValue() : ""; }

  T& target;
};

template <typename T>
struct LowerBoundArgumentOption : ArgumentOption<T> {
  LowerBoundArgumentOption(T& target_, char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_, std::string argName, T minValue_)
      : ArgumentOption<T>(target_, shortName_, std::move(longName_), std::move(desc), std::move(group_), isAdvanced_, std::move(argName)), minValue(minValue_) {}

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
  LowerUpperBoundArgumentOption(T& target_, char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_, std::string argName, T minValue_, T maxValue_)
      : LowerBoundArgumentOption<T>(target_, shortName_, std::move(longName_), std::move(desc), std::move(group_), isAdvanced_, std::move(argName), minValue_), maxValue(maxValue_) {}

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
  ArgumentOption(bool& target_, char shortName_, std::string longName_, std::string desc, std::string group_, bool isAdvanced_)
      : Option(shortName_, std::move(longName_), std::move(desc), std::move(group_), isAdvanced_, false), target(target_) {}

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
  std::string printValue() const override { return io::stringify(target); }
  std::string printNumericValue() const override { return ""; }

  bool& target;
};

struct TargetBase {
  TargetBase(std::string variableName_, std::string desc, std::string group_, bool isAdvanced_)
      : variableName(std::move(variableName_)), description(std::move(desc)), group(std::move(group_)), isOptionalVector(false), isAdvanced(isAdvanced_) {}

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
  Target(T& target_, std::string variableName_, std::string desc, std::string group_, bool isAdvanced_)
      : TargetBase(std::move(variableName_), std::move(desc), std::move(group_), isAdvanced_), target(target_) {}

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

  void setJsonParameters(std::string jsonParameters) { m_jsonParameters = std::move(jsonParameters); }
  void setJsonParametersProvider(std::function<std::string()> jsonParametersProvider) { m_jsonParametersProvider = std::move(jsonParametersProvider); }

  std::vector<ParsedOption> getUsedOptionArguments() const;

  unsigned int numRequiredArguments() const { return m_nonOptionArguments.size() - m_numOptionalNonOptionArguments; }

private:
  void exitWithUsage(bool showAdvanced) const;
  void exitWithVersionInformation() const;
  void exitWithError(const std::string& message) const;
  void exitWithCompletion() const;
  void exitWithJsonParameters() const;
  void printBashCompletion() const;
  void printZshCompletion() const;

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
  std::string m_completionShell;
  std::string m_jsonParameters;
  std::function<std::string()> m_jsonParametersProvider;

  unsigned int m_numOptionalNonOptionArguments = 0;
};

} // namespace infomap

#endif // PROGRAM_INTERFACE_H_

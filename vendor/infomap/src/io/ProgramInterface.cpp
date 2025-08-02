/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "ProgramInterface.h"
#include "../utils/Log.h"

#include <iostream>
#include <cstdlib>
#include <map>
#include <utility>

namespace infomap {

const std::string ArgType::integer = "integer";
const std::string ArgType::number = "number";
const std::string ArgType::string = "string";
const std::string ArgType::path = "path";
const std::string ArgType::probability = "probability";
const std::string ArgType::option = "option";
const std::string ArgType::list = "list";

const std::unordered_map<std::string, char> ArgType::toShort = {
  { "integer", 'n' },
  { "number", 'f' },
  { "string", 's' },
  { "path", 'p' },
  { "probability", 'P' },
  { "option", 'o' },
  { "list", 'l' },
};

ProgramInterface::ProgramInterface(std::string name, std::string shortDescription, std::string version)
    : m_programName(std::move(name)),
      m_shortProgramDescription(std::move(shortDescription)),
      m_programVersion(std::move(version))
{
  addIncrementalOptionArgument(m_displayHelp, 'h', "help", "Prints this help message. Use -hh to show advanced options.", "About");
  addOptionArgument(m_displayVersion, 'V', "version", "Display program version information.", "About");
  addOptionArgument(m_printJsonParameters, "print-json-parameters", "Print Infomap parameters in JSON.", "About").setHidden(true);
}

void ProgramInterface::exitWithUsage(bool showAdvanced) const
{
  Log() << "Name:\n";
  Log() << "        " << m_programName << " - " << m_shortProgramDescription << '\n';
  Log() << "\nUsage:\n";
  Log() << "        " << m_executableName;
  for (auto& nonOptionArgument : m_nonOptionArguments)
    if (showAdvanced || !nonOptionArgument->isAdvanced)
      Log() << " " << nonOptionArgument->variableName;
  if (!m_optionArguments.empty())
    Log() << " [options]";
  Log() << '\n';

  if (!m_programDescription.empty())
    Log() << "\nDescription:\n        " << m_programDescription << '\n';

  for (auto& nonOptionArgument : m_nonOptionArguments)
    if (showAdvanced || !nonOptionArgument->isAdvanced)
      Log() << "\n[" << nonOptionArgument->variableName << "]\n    " << nonOptionArgument->description << '\n';

  if (!m_optionArguments.empty())
    Log() << "\n[options]\n";

  // First stringify the options part to get the maximum length
  std::deque<std::string> optionStrings(m_optionArguments.size());
  std::string::size_type maxLength = 0;
  for (unsigned int i = 0; i < m_optionArguments.size(); ++i) {
    auto& opt = *m_optionArguments[i];
    bool haveShort = opt.shortName != '\0';
    std::string optArgShort = opt.requireArgument ? (io::Str() << "<" << ArgType::toShort.at(opt.argumentName) << ">") : opt.incrementalArgument ? "[+]"
                                                                                                                                                 : std::string(3, ' ');
    std::string optArgLong = opt.requireArgument ? (io::Str() << "<" << opt.argumentName << ">") : opt.incrementalArgument ? "[+]"
                                                                                                                           : std::string(3, ' ');
    std::string shortOption = haveShort ? (io::Str() << "  -" << opt.shortName << optArgShort) : std::string(7, ' ');
    optionStrings[i] = io::Str() << shortOption << " --" << opt.longName << " " << optArgLong;
    if (optionStrings[i].length() > maxLength)
      maxLength = optionStrings[i].length();
  }

  std::vector<std::string> groups { "About" };
  for (auto& group : m_groups) {
    if (group != "About")
      groups.push_back(group);
  }
  if (m_groups.empty())
    groups.emplace_back("All");

  for (const auto& group : groups) {
    if (group != "All") {
      Log() << "\n"
            << group << "\n";
      Log() << std::string(group.length(), '-') << "\n";
    }
    for (unsigned int i = 0; i < m_optionArguments.size(); ++i) {
      auto& opt = *m_optionArguments[i];
      if (group == "All" || opt.group == group) {
        std::string::size_type numSpaces = maxLength + 3 - optionStrings[i].length();
        if (showAdvanced || !opt.isAdvanced) {
          Log() << optionStrings[i] << std::string(numSpaces, ' ') << opt.description;
          if (!opt.printNumericValue().empty())
            Log() << " (Default: " << opt.printNumericValue() << ")";
          Log() << "\n";
        }
      }
    }
  }
  Log() << '\n';
  std::exit(0);
}

void ProgramInterface::exitWithVersionInformation() const
{
  Log() << m_programName << " version " << m_programVersion;
#ifdef _OPENMP
  Log() << " compiled with OpenMP";
#endif
  Log() << '\n';
  Log() << "See www.mapequation.org for terms of use.\n";
  std::exit(0);
}

void ProgramInterface::exitWithError(const std::string& message) const
{
  Log() << m_programName << " version " << m_programVersion;
#ifdef _OPENMP
  Log() << " compiled with OpenMP";
#endif
  Log() << std::endl;
  std::cerr << message << std::endl;
  Log() << "Usage: " << m_executableName;
  for (auto& nonOptionArgument : m_nonOptionArguments)
    if (!nonOptionArgument->isAdvanced)
      Log() << " " << nonOptionArgument->variableName;
  if (!m_optionArguments.empty())
    Log() << " [options]";
  Log() << ". Run with option '-h' for more information.\n";
  std::exit(1);
}

std::string toJson(const std::string& key, const std::string& value)
{
  return io::Str() << '"' << key << "\": \"" << value << '"';
}

std::string toJson(const std::string& key, bool value)
{
  return io::Str() << '"' << key << "\": " << (value ? "true" : "false");
}

template <typename Value>
std::string toJson(const std::string& key, Value value)
{
  return io::Str() << '"' << key << "\": " << value;
}

std::string toJson(const Option& opt)
{
  return io::Str() << "{ "
                   << toJson("long", std::string(io::Str() << "--" << opt.longName)) << ", "
                   << toJson("short", opt.shortName != '\0' ? std::string(io::Str() << "-" << opt.shortName) : "") << ", "
                   << toJson("description", opt.description) << ", "
                   << toJson("group", opt.group) << ", "
                   << toJson("required", opt.requireArgument) << ", "
                   << toJson("advanced", opt.isAdvanced) << ", "
                   << toJson("incremental", opt.incrementalArgument) << ", "
                   << (opt.requireArgument
                           ? (io::Str() << toJson("longType", opt.argumentName) << ", "
                                        << toJson("shortType", std::string(1, ArgType::toShort.at(opt.argumentName))) << ", "
                                        << toJson("default", opt.printValue()))
                           : toJson("default", false))
                   << " }";
}

void ProgramInterface::exitWithJsonParameters() const
{
  Log() << "{\n  \"parameters\": [\n";

  for (unsigned int i = 0; i < m_optionArguments.size(); ++i) {
    auto& opt = *m_optionArguments[i];
    if (opt.hidden)
      continue;
    Log() << "    " << toJson(opt);
    if (i < m_optionArguments.size() - 1) {
      Log() << ",\n";
    } else {
      Log() << "\n";
    }
  }
  Log() << "  ]\n}";

  std::exit(0);
}

void ProgramInterface::parseArgs(const std::string& args)
{
  // Map the options on short and long name, and check for duplication
  std::map<char, Option*> shortOptionMap;
  std::map<std::string, Option*> longOptionMap;
  for (auto& optionArgument : m_optionArguments) {
    auto& opt = *optionArgument;
    if (opt.shortName != '\0') {
      auto it = shortOptionMap.find(opt.shortName);
      if (it != shortOptionMap.end())
        throw std::runtime_error(io::Str() << "Duplication of option '" << opt.shortName << "'");
      shortOptionMap.insert(std::make_pair(opt.shortName, &opt));
    }

    auto it = longOptionMap.find(opt.longName);
    if (it != longOptionMap.end())
      throw std::runtime_error(io::Str() << "Duplication of option \"" << opt.longName << "\"");
    longOptionMap.insert(std::make_pair(opt.longName, &opt));
  }

  // Split the flags on whitespace
  std::vector<std::string> flags;
  std::istringstream argStream(args);

  {
    std::string arg;
    while (!(argStream >> arg).fail())
      flags.push_back(arg);
  }

  std::deque<std::string> nonOpts;
  try {
    for (unsigned int i = 0; i < flags.size(); ++i) {
      bool flagValue = true;
      unsigned int numArgsLeft = flags.size() - i - 1;

      const std::string& arg = flags[i];
      if (arg.length() == 0)
        throw std::runtime_error("Illegal argument ''");

      if (arg[0] != '-') {
        nonOpts.push_back(arg);
      } else {
        if (arg.length() < 2)
          throw std::runtime_error("Illegal argument '-'");

        if (arg[1] == '-') {
          // Long option
          if (arg.length() < 3)
            throw std::runtime_error("Illegal argument '--'");
          std::string longOpt = arg.substr(2);
          auto it = longOptionMap.find(longOpt);
          if (it == longOptionMap.end()) {
            // Unrecognized option, check if it negates a recognised option with the '--no-' prefix
            if (longOpt.compare(0, 3, "no-") == 0 && longOptionMap.find(std::string(longOpt, 3)) != longOptionMap.end()) {
              longOpt = std::string(longOpt, 3);
              it = longOptionMap.find(longOpt);
              flagValue = false;
            } else {
              throw std::runtime_error(io::Str() << "Unrecognized option: '--" << longOpt << "'");
            }
          }
          auto& opt = *it->second;
          if (!opt.requireArgument || opt.incrementalArgument)
            opt.set(flagValue);
          else {
            if (numArgsLeft == 0)
              throw std::runtime_error(io::Str() << "Option '" << opt.longName << "' requires argument");
            ++i;
            if (!opt.parse(flags[i]))
              throw std::runtime_error(io::Str() << "Cannot parse '" << flags[i] << "' as argument to option '" << opt.longName << "'. ");
          }
        } else {
          // Short option(s)
          for (unsigned int j = 1; j < arg.length(); ++j) {
            char o = arg[j];
            unsigned int numCharsLeft = arg.length() - j - 1;
            auto it = shortOptionMap.find(o);
            if (it == shortOptionMap.end())
              throw std::runtime_error(io::Str() << "Unrecognized option: '-" << o << "'");
            auto& opt = *it->second;
            if (!opt.requireArgument || opt.incrementalArgument)
              opt.set(flagValue);
            else {
              std::string optArg;
              if (numCharsLeft > 0) {
                optArg = arg.substr(j + 1);
                j = arg.length() - 1;
              } else if (numArgsLeft) {
                ++i;
                optArg = flags[i];
              } else
                throw std::runtime_error(io::Str() << "Option '" << opt.longName << "' requires argument");

              if (!opt.parse(optArg))
                throw std::runtime_error(io::Str() << "Cannot parse '" << optArg << "' as argument to option '" << opt.longName << "'. ");
            }
          }
        }
      }
      if (m_displayHelp > 0)
        exitWithUsage(m_displayHelp > 1);
      if (m_displayVersion)
        exitWithVersionInformation();
      if (m_printJsonParameters)
        exitWithJsonParameters();
    }
  } catch (std::exception& e) {
    exitWithError(e.what());
  }

  if (nonOpts.size() < numRequiredArguments())
    exitWithError("Missing required arguments.");

  unsigned int i = 0;
  unsigned int numVectorArguments = nonOpts.size() - (m_nonOptionArguments.size() - 1);
  while (!nonOpts.empty()) {
    std::string arg = nonOpts.front();
    nonOpts.pop_front();
    if (m_nonOptionArguments[i]->isOptionalVector && numVectorArguments == 0)
      ++i;
    if (!m_nonOptionArguments[i]->parse(arg))
      exitWithError("Argument error.");
    if (!m_nonOptionArguments[i]->isOptionalVector || --numVectorArguments == 0)
      ++i;
  }
}

std::vector<ParsedOption> ProgramInterface::getUsedOptionArguments() const
{
  std::vector<ParsedOption> opts;
  unsigned int numFlags = m_optionArguments.size();
  for (unsigned int i = 0; i < numFlags; ++i) {
    auto& opt = *m_optionArguments[i];
    if (opt.used && opt.longName != "negate-next")
      opts.emplace_back(opt);
  }
  return opts;
}

} // namespace infomap

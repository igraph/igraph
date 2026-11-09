#!/usr/bin/env Rscript
# Generate interfaces/functions-categories.yaml from doc/*.xxml files
# Usage: Rscript tools/generate_function_categories.R

library(yaml)

# Paths
xxml_dir <- "doc"
functions_yaml <- "interfaces/functions.yaml"
output_yaml <- "interfaces/functions-categories.yaml"

# Read function names from functions.yaml
functions <- yaml.load_file(functions_yaml)
function_names <- names(functions)
cat("Total functions in functions.yaml:", length(function_names), "\n")

# Parse all .xxml files
xxml_files <- list.files(xxml_dir, pattern = "\\.xxml$", full.names = TRUE)
cat("Found", length(xxml_files), ".xxml files\n")

# Build mapping: function_name -> list(category, subcategory)
func_to_info <- list()

for (f in xxml_files) {
  category <- tools::file_path_sans_ext(basename(f))
  content <- readLines(f, warn = FALSE)

  # Track a stack of section ids to find the innermost one
  section_id_stack <- character(0)

  for (line in content) {
    # Opening <section ...> — push id (or NA) onto stack
    if (grepl("<section", line)) {
      id_match <- regmatches(line, regexpr('(?<=id=")[^"]+(?=")', line, perl = TRUE))
      section_id_stack <- c(section_id_stack, if (length(id_match) == 1) id_match else NA_character_)
    }

    # Closing </section> — pop from stack
    if (grepl("</section>", line)) {
      if (length(section_id_stack) > 0)
        section_id_stack <- section_id_stack[-length(section_id_stack)]
    }

    # doxrox-include — record function with current subcategory
    name_match <- regmatches(line, regexpr("(?<=<!-- doxrox-include )\\S+(?= -->)", line, perl = TRUE))
    if (length(name_match) == 1 && name_match %in% function_names) {
      name <- name_match
      if (!is.null(func_to_info[[name]])) {
        warning("Function ", name, " appears in multiple categories: ",
                func_to_info[[name]]$category, " and ", category)
      }
      # Use innermost section id that is not NA
      valid_ids <- section_id_stack[!is.na(section_id_stack)]
      subcategory <- if (length(valid_ids) > 0) valid_ids[length(valid_ids)] else NA_character_
      func_to_info[[name]] <- list(category = category, subcategory = subcategory)
    }
  }
}

cat("Functions with a category:", length(func_to_info), "\n")

# Undocumented functions
undocumented <- setdiff(function_names, names(func_to_info))
cat("Undocumented functions:", length(undocumented), "\n")

# Build output data frame: sort by category, subcategory, name
categorized_df <- data.frame(
  name        = names(func_to_info),
  category    = sapply(func_to_info, `[[`, "category"),
  subcategory = sapply(func_to_info, `[[`, "subcategory"),
  stringsAsFactors = FALSE
)
categorized_df <- categorized_df[
  order(categorized_df$category, categorized_df$subcategory, categorized_df$name),
]

# Write YAML output
lines <- c("# Functions ordered by category")

for (i in seq_len(nrow(categorized_df))) {
  entry <- c(paste0(categorized_df$name[i], ":"),
             paste0("    category: ", categorized_df$category[i]))
  if (!is.na(categorized_df$subcategory[i])) {
    entry <- c(entry, paste0("    subcategory: ", categorized_df$subcategory[i]))
  }
  lines <- c(lines, entry, "")
}

# Append undocumented entries
if (length(undocumented) > 0) {
  lines <- c(lines, "# Functions not found in any .xxml file")
  for (name in sort(undocumented)) {
    lines <- c(lines,
               paste0(name, ":"),
               "    category: undocumented",
               "")
  }
}

writeLines(lines, output_yaml)
cat("Written to", output_yaml, "\n")

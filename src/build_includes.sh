#!/bin/bash

# Define the top-level include file
top_include="src/include.hpp"

# Start with #pragma once in the top-level file
printf "#pragma once\n" > "$top_include"

# Special case: Dirac Operator:
search_dir="src/DiracOperator/Operators"
output_file="src/DiracOperator/Operators.hpp"

# Create or overwrite the output file with the pragma directive
printf "#pragma once\n" > "$output_file"
# Loop through all .hpp files in the directory
for header in "$search_dir"/*.hpp; do
    # Get the relative path from 'src/' to the file
    rel_header="${header#src/}"

    # Add an #include for each header file with the full relative path
    printf "#include \"%s\"\n" "$rel_header" >> "$output_file"
done

# Loop through all directories inside src/
for d in src/*/ ; do
    # Remove only "src/" but keep the subdirectory name
    rel_dir="${d#src/}"          # e.g., "subdir1/"
    rel_dir="${rel_dir%/}"       # Remove trailing slash -> "subdir1"
    
    t_filename="${d}include.hpp"  # Full path for the subdir include file
    rel_t_filename="${rel_dir}/include.hpp"  # Relative path for use in #include

    # Check if the directory contains any .hpp files
    if compgen -G "${d}*.hpp" > /dev/null; then
        # Create or overwrite include.hpp with the pragma directive
        printf "#pragma once\n" > "$t_filename"

        # Loop through all .hpp files in the directory
        for header in "${d}"*.hpp; do
            filename=$(basename "$header")

            # Skip "include.hpp"
            if [[ "$filename" != "include.hpp" ]]; then
                printf "#include \"%s/%s\"\n" "$rel_dir" "$filename" >> "$t_filename"
            fi
        done

        # Add this include file to the top-level include file
        printf "#include \"%s\"\n" "$rel_t_filename" >> "$top_include"
    fi
done
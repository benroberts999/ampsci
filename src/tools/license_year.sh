#!/bin/bash
# Updates the copyright year in LICENSE and copies to doc/LICENSE.md.

year=$(date +%Y)
sed -i.bak -E "s/(Copyright \(c\) )([0-9]{4})(-[0-9]{4})?/\1\2-${year}/" LICENSE
rm -f LICENSE.bak
cp ./LICENSE ./doc/LICENSE.md
echo "Updated license year to ${year}"

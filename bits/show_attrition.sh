#!/bin/bash
echo "First 10 rows of group count files for major processing stages"
for f in $((find . -name '*.groupstats' -exec ls -1rt "{}" +;) | grep OUTPUT); do echo $f; head -n 10 $f; done


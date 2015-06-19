#!/usr/bin/env python
## Otu000390       18      Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Blautia(100);
import sys
for i_line,line in enumerate(sys.stdin):
    if i_line>0 and len(line.split(";")) == 7:
        line = line.strip()+"unclassified(100);\n"
    sys.stdout.write(line)


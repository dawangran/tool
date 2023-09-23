
###################################################
# tilte: bam change tag for bgi
# date:20230922
# author: dawn
###################################################

import sys
import pysam
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
outsam = pysam.AlignmentFile(sys.argv[2], "wb",header=dict(samfile.header))

for i in samfile:
    try:
        if i.has_tag('UB'):
            RCB = i.get_tag('CB')
            RDB = i.get_tag('DB')
            i.set_tag('DB',RCB)
            i.set_tag('CB',RDB)
            outsam.write(i)
    except KeyError:
        continue

outsam.close()
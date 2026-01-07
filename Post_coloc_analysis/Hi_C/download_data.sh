#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

#wget ftp://ftp.broadinstitute.org/outgoing/lincRNA/average_hic/average_hic.v2.191020.tar.gz

# The folder I am sharing contains five blood cell lines including K562 and a control keratinocyte cell line.
# However, it is all mapped to hg19. Each folder contains three different sets of interactions called at different FDR cutoffs.
# I recommend you start with the most lenient cutoff.

# I couldn't download via the link using wget.


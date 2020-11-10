import sys
import os
from shutil import copyfile

directory = (sys.argv[1])
genome = (sys.argv[2])


species = directory.replace("/", "")


hub_text = f"hub {species}\n" \
           f"shortLabel {species} Popgen Stats\n" \
           f"longLabel {species} aligned to {genome} Popgen Stats\n" \
           f"genomesFile genomes.txt\n" \
           f"email cmirchan@ucsc.edu"

hub_file = os.path.join(directory, "hub.txt")

with open(hub_file, 'w') as h:
    h.write(hub_text)


genome_text = f"genome {genome}\n" \
              f"trackDb {genome}/trackDb.txt"

genome_file = os.path.join(directory, "genomes.txt")

with open(genome_file, 'w') as g:
    g.write(genome_text)

trackdb = os.path.join(os.getcwd(), "helper_files", "trackDb.txt")

if os.path.exists(os.path.join(directory, genome)):
    copyfile(trackdb, os.path.join(directory, genome,"trackDb.txt"))
else:
    os.mkdir(os.path.join(directory, genome))
    copyfile(trackdb, os.path.join(directory, genome,"trackDb.txt"))
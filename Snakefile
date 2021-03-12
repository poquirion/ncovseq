from packaging import version
from snakemake.logging import logger
import os, sys
from shutil import which
import datetime

#
# Verify that required versions of dependencies are installed.
#
MIN_AUGUR_VERSION = "7.0.2"

try:
    from augur.__version__ import __version__ as augur_version
except ModuleNotFoundError:
    logger.error("ERROR: Could not find augur. Follow installation instructions at https://nextstrain.org/docs/ and try again.")
    sys.exit(1)

if version.parse(augur_version) < version.parse(MIN_AUGUR_VERSION):
    logger.error("ERROR: Found version '%s' of augur, but version '%s' or greater is required" % (augur_version, MIN_AUGUR_VERSION))
    sys.exit(1)

SHELL_COMMANDS_NEEDED = ["augur", "iqtree", "mafft"]
for sh_cmd in SHELL_COMMANDS_NEEDED:
    if not which(sh_cmd):
        logger.error(f"ERROR: `{sh_cmd}` is not available as a shell command. Please follow installation instructions at https://nextstrain.org/docs/ and try again.")
        sys.exit(1)


def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

# For information on how to run 'regions' runs, see Snakefile_Regions

# Add new regions here!
REGIONS = ["_quebec","_global"]

TRAIT_REZ = ('location','province')

wildcard_constraints:
    region = "|".join(REGIONS) + "||",
    date = "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

configfile: "config/Snakefile.yaml"


# simple rule to call snakemake for outsider users
rule all:
    input:
        auspice_json = "auspice/ncov.json",
        tip_frequencies_json = "auspice/ncov_tip-frequencies.json",


include: "workflow/snakemake_rules/organize_meta.smk"
include: "workflow/snakemake_rules/select_align_fasta.smk"
include: "workflow/snakemake_rules/build_tree.smk"
include: "workflow/snakemake_rules/prepare_visu.smk"


rule tree_only:
    input:
        tree = rules.refine.output.tree,
        alignement = rules.aggregate_alignments.output.alignment



rule redo_maps:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        """
        rm -f {rules.fix_colorings.output.auspice_json}
        rm -f {rules.tip_frequencies.output.tip_frequencies_json}
        rm -f {rules.export.output.auspice_json}
        rm -f {rules.lat_long.output.lat_long}
        """


rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "find {params} -type f -not -name 'nextstrain_exclude*' -print0 | xargs --null -P 33 -L 1  rm  "



rule reprod_zip:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

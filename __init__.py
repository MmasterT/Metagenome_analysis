import pkg_resources

__title__ = "yaragenome"
__author__ = "Mariano Olivera Fedi"
__email__ = "mariano.Olivera-Fedi@earlham.ac.uk"
__copyright__ = "Copyright 2022 Earlham Institute"
__version__ = pkg_resources.require("yaragenome")[0].version

DEFAULT_CONFIG_FILE = pkg_resources.resource_filename("yaragenome.etc", "run_config.yaml")
DEFAULT_HPC_CONFIG_FILE = pkg_resources.resource_filename(
    "yaragenome.etc", "hpc_config.json"
)
DEFAULT_SAMPLE_CSV_FILE = pkg_resources.resource_filename("yaragenome.etc", "samples.csv")

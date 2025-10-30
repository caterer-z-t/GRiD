# grid/utils/helper_dir/load_flags_from_yaml.py
# In[1]: Imports
import yaml
from pathlib import Path
from typing import Set

# In[2]: Function to load read count flags from config file
def load_flags(config_file: str, parameter: str) -> Set[int]:
    """
    Load read count flags from a YAML configuration file.

    Args:
        config_file (str): Path to the YAML configuration file.
        parameter (str): The parameter name to look for in the YAML file.

    Returns:
        Set[int]: A set of read count flags.
    """

    config_file = Path(config_file).expanduser()
    with open(config_file, "r") as f:
        cfg = yaml.safe_load(f)

    # Support both possible formats
    if parameter in cfg:
        flags = cfg[parameter]['flags']
    else:
        raise ValueError(f"No '{parameter}' or 'read-count.{parameter}' found in {config_file}")

    return set(int(f) for f in flags)
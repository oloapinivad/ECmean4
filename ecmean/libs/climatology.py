"""
Common functions for reference and climatology creation.
"""

import logging

loggy = logging.getLogger(__name__)


# Climatology file prefixes
CLIMATOLOGY_PREFIXES = ["climate_variance_", "climate_average_", "seasons_variance_", "seasons_average_"]

# supported climatology
SUPPORTED_REFERENCE = ["EC23", "EC26-PDAY", "EC26-CMIP", "EC26-HIST"]
SUPPORTED_CLIMATOLOGY = ["EC23", "EC24", "EC26-HIST", "EC26-CMIP"]

"""GWAS (Genome-Wide Association Studies) module.

This module provides classes for performing GWAS analysis using different models:
- Preparatory steps (LD pruning, PCA decomposition)
- Fixed effects model (Generalized Linear Model)
- Random effects model (Generalized Linear Mixed Model)
"""

from ideal_genom.gwas.preparatory import Preparatory
from ideal_genom.gwas.gen_linear_model import GWAS_GLM
from ideal_genom.gwas.gen_linear_mix_model import GWAS_GLMM

__all__ = [
    'Preparatory',
    'GWAS_GLM',
    'GWAS_GLMM',
]

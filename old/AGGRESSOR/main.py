#!/usr/bin/env python3
"""
AGGRESSOR: Aggregation-Guided Generation of REgion-Specific Substitution ORiented mutations

Main entry point for the AGGRESSOR pipeline.
Delegates all logic to the CLI module after ensuring clean logger state.

Usage:
    python main.py <input_file> [options]
    python -m aggressor <input_file> [options]

References:
- Rousseau et al., J Mol Biol 2006 (gatekeeper hypothesis)
- Beerten et al., FEBS Lett 2012 (APR boundary effects)
- Tartaglia et al., J Mol Biol 2008 (aggregation propensity scale)
"""
import logging

from CLI import main

if __name__ == '__main__':
    # Prevent duplicate log handlers when running as script
    # This is especially important in development environments
    # where auto-reload might re-execute the module
    logging.getLogger().handlers.clear()

    main()
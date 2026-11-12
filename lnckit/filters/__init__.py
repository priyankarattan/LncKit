"""Transcript filters that progressively refine candidates into lncRNAs.

Each filter takes a set of transcript IDs (plus whatever evidence file it needs)
and returns the surviving subset, so filters compose cleanly and every step's
retained count is recorded for the report.
"""

"""Stage orchestrators for the upstream (read-level) half of the pipeline.

Each function takes the RunContext, runs the configured backend for that stage
across all samples, records a checkpoint, and returns the paths the next stage
needs. They deliberately contain no tool-specific command logic - that lives in
the wrappers - so swapping a backend never touches these modules.
"""

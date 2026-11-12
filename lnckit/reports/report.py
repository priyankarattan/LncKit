"""Render run_metadata.json into an HTML report (and optionally PDF)."""
from __future__ import annotations
import json
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.report")
_TPL_DIR = Path(__file__).parent / "templates"


class _M:  # lightweight attribute view over the metadata dict
    def __init__(self, d): self.__dict__.update(d)


def build_report(output_dir: Path, to_pdf: bool = False) -> Path:
    output_dir = Path(output_dir)
    meta_path = output_dir / "reports" / "run_metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"No run metadata at {meta_path}; run the pipeline first.")
    meta = json.loads(meta_path.read_text())

    env = Environment(loader=FileSystemLoader(str(_TPL_DIR)),
                      autoescape=select_autoescape(["html"]))
    tpl = env.get_template("report.html")
    steps = [_M(s) for s in meta.get("steps", [])]
    m = _M({**meta, "steps": steps, "final_counts": meta.get("final_counts", {})})
    html = tpl.render(m=m, config_yaml=yaml.safe_dump(meta.get("config", {}), sort_keys=False))

    html_out = output_dir / "reports" / "report.html"
    html_out.write_text(html)
    log.info("HTML report: %s", html_out)

    if to_pdf:
        try:
            from weasyprint import HTML
            pdf_out = output_dir / "reports" / "report.pdf"
            HTML(string=html).write_pdf(str(pdf_out))
            log.info("PDF report: %s", pdf_out)
        except Exception as e:  # pragma: no cover
            log.warning("PDF rendering skipped (%s). Install lnckit[reports].", e)
    return html_out

#!/usr/bin/env python

import importlib.util
import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
PACKAGE_DIR = ROOT / "output" / "frontiers_submission_package"
SUBMISSION_DIR = WORKSPACE / "submission_frontiers_oncology_revision"
BUILDER = ROOT / "scripts" / "15_build_frontiers_submission_package.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


builder = load_module("frontiers_builder", BUILDER)


def pass_fail(label, ok, detail):
    status = "PASS" if ok else "FAIL"
    return "[%s] %s: %s" % (status, label, detail)


def main():
    lines = []
    ok_all = True

    submission_tables_ok = (
        (SUBMISSION_DIR / "tables" / "Table2_external_validation_summary.tsv").exists()
    )
    lines.append(pass_fail("Submission targets present", submission_tables_ok, "submission table files used as the final comparison target"))
    ok_all &= submission_tables_ok

    readme_text = (PACKAGE_DIR / "README_submission_package.md").read_text(encoding="utf-8")
    readme_ok = ("six independent external validation cohorts" not in readme_text and "Choueiri2016" not in readme_text and "GSE29609" in readme_text)
    lines.append(pass_fail("Package README", readme_ok, "README uses the revised Frontiers framing"))
    ok_all &= readme_ok

    table2 = pd.read_csv(PACKAGE_DIR / "tables" / "Table2_external_validation_summary.tsv", sep="\t", dtype=str)
    expected_table2 = pd.DataFrame(builder.TABLE2_ROWS)[table2.columns.tolist()].astype(str)
    table2_ok = table2.equals(expected_table2)
    lines.append(pass_fail("Table 2 retained validation cohorts", table2_ok, "expected 4 retained cohorts, observed %d" % table2.shape[0]))
    ok_all &= table2_ok

    table4 = pd.read_csv(PACKAGE_DIR / "tables" / "Table4_orientation_sensitivity.tsv", sep="\t", dtype=str)
    expected_table4 = pd.DataFrame(builder.TABLE4_ROWS).astype(str)
    table4_ok = table4.equals(expected_table4)
    lines.append(pass_fail("Table 4 orientation rows", table4_ok, "expected 8 rows, observed %d" % table4.shape[0]))
    ok_all &= table4_ok

    table5 = pd.read_csv(PACKAGE_DIR / "tables" / "Table5_clinicopathologic_benchmarking.tsv", sep="\t", dtype=str)
    expected_table5 = pd.DataFrame(builder.TABLE5_EXPECTED).astype(str)
    table5_ok = table5.equals(expected_table5)
    lines.append(pass_fail("Table 5 benchmarking cohorts", table5_ok, "expected TCGA-KIRC plus exploratory E-MTAB-1980 only"))
    ok_all &= table5_ok

    table3 = pd.read_csv(PACKAGE_DIR / "tables" / "Table3_final_33_gene_signature.tsv", sep="\t", dtype=str)
    coef_ok = table3["Coefficient"].str.fullmatch(r"-?\d+\.\d{6}").all() and table3["Absolute coefficient"].str.fullmatch(r"\d+\.\d{6}").all()
    lines.append(pass_fail("Table 3 coefficient formatting", coef_ok, "coefficients formatted to exactly 6 decimal places"))
    ok_all &= coef_ok

    figure_names = {
        "Figure1_strobe_flow_diagram.tiff",
        "Figure2_study_workflow.tiff",
        "Figure3_external_kaplan_meier.tiff",
        "Figure4_discrimination_summary.tiff",
        "Figure5_robustness_diagnostics.tiff",
        "Figure6_signature_characteristics.tiff",
    }
    observed_names = {path.name for path in (PACKAGE_DIR / "figures").glob("*.tiff")}
    figures_ok = figure_names.issubset(observed_names)
    lines.append(pass_fail("Figure TIFF set", figures_ok, "expected 6 manuscript figure TIFF files"))
    ok_all &= figures_ok

    report_path = PACKAGE_DIR / "consistency_check_frontiers.txt"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    if not ok_all:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

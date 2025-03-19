#!/usr/bin/env python3
import argparse
import sqlite3
import json
import plotly.graph_objects as go
import numpy as np


def parse_float_or_nan(x):
    """Convert string values to float or NaN."""
    if isinstance(x, str) and x.lower() == "nan":
        return float("nan")
    return float(x)


def query_allele_data(db_path, repeat_id):
    """
    Queries the SQLite database for allele data using repeat_id.
    """
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("SELECT data_json FROM locus_data WHERE repeat_id = ?", (repeat_id,))
    result = cur.fetchone()
    conn.close()

    if not result:
        print(f"[WARNING] No data found for repeat_id: {repeat_id}")
        return None, None, None

    data = json.loads(result[0])

    # Extract relevant fields
    dosage_dict = json.loads(data.get("sample_count_per_summed_length", "{}"))
    mean_col_name = next((c for c in data.keys() if c.startswith("mean_")), None)
    mean_dict = json.loads(data.get(mean_col_name, "{}")) if mean_col_name else {}
    ci_dict = json.loads(data.get("summed_length_0.05_alpha_CI", "{}"))

    return dosage_dict, mean_dict, ci_dict


def filter_allele_data(dosage_dict, mean_dict, ci_dict, count_threshold, max_ci_range, max_relative_ci_range):
    """
    Applies threshold-based filtering to allele data.
    """
    filtered_dosage = {}
    filtered_mean = {}
    filtered_ci = {}

    for allele in dosage_dict.keys():
        count_val = parse_float_or_nan(dosage_dict[allele])
        if count_val < count_threshold:
            continue

        lower = parse_float_or_nan(ci_dict.get(allele, [None, None])[0])
        upper = parse_float_or_nan(ci_dict.get(allele, [None, None])[1])
        mean_val = parse_float_or_nan(mean_dict.get(allele, None))

        if np.isnan(lower) or np.isnan(upper) or np.isnan(mean_val):
            continue

        # Apply CI filtering
        if max_ci_range is not None and (upper - lower) > max_ci_range:
            continue
        if max_relative_ci_range is not None and ((upper - lower) / mean_val) > max_relative_ci_range:
            continue

        filtered_dosage[allele] = dosage_dict[allele]
        filtered_mean[allele] = mean_dict[allele]
        filtered_ci[allele] = ci_dict[allele]

    return filtered_dosage, filtered_mean, filtered_ci


def generate_figure_plotly(dosage_dict, mean_dict, ci_dict):
    """
    Creates a Plotly figure using the processed allele data.
    Defaults to black & white (B/W) with error bars.
    """
    sorted_alleles = sorted(dosage_dict.keys(), key=float)
    if not sorted_alleles:
        print("[ERROR] No valid allele data found.")
        return None

    ci_lower = [ci_dict.get(str(a), [None, None])[0] for a in sorted_alleles]
    ci_upper = [ci_dict.get(str(a), [None, None])[1] for a in sorted_alleles]
    mean_vals = [mean_dict.get(str(a), None) for a in sorted_alleles]

    fig = go.Figure()

    # ✅ Default to Black & White (B/W) error bars
    error_y = [ci_upper[i] - mean_vals[i] for i in range(len(sorted_alleles))]
    error_y_minus = [mean_vals[i] - ci_lower[i] for i in range(len(sorted_alleles))]
    fig.add_trace(go.Scatter(
        x=sorted_alleles, y=mean_vals, mode="lines+markers",
        error_y=dict(type="data", array=error_y, arrayminus=error_y_minus, visible=True),
        line=dict(color="black", width=3), marker=dict(color="black", size=8), name="95% CI"
    ))

    fig.update_layout(
        xaxis_title="Sum of allele lengths (repeat copies)",
        yaxis_title="Phenotype Value",
        showlegend=True
    )

    return fig



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-path", required=True, help="Path to the SQLite database file.")
    parser.add_argument("--repeat-id", required=True, help="Repeat ID to query.")
    parser.add_argument("--output-dir", required=True, help="Directory to save the output.")
    parser.add_argument("--count-threshold", type=float, default=100, help="Minimum allele count threshold.")
    parser.add_argument("--max-ci-range", type=float, default=None, help="Max absolute CI range.")
    parser.add_argument("--max-relative-ci-range", type=float, default=None, help="Max relative CI range.")

    args = parser.parse_args()

    dosage_dict, mean_dict, ci_dict = query_allele_data(args.db_path, args.repeat_id)
    if not dosage_dict:
        print("[ERROR] No valid data found. Exiting.")
        return

    dosage_dict, mean_dict, ci_dict = filter_allele_data(
        dosage_dict, mean_dict, ci_dict, args.count_threshold, args.max_ci_range, args.max_relative_ci_range
    )

    fig = generate_figure_plotly(dosage_dict, mean_dict, ci_dict)

    if fig:
        output_path = f"{args.output_dir}/{args.repeat_id}.html"
        fig.write_html(output_path)
        fig.show()  # Display the plot as a pop-up
        print(f"Saved interactive plot to {output_path}")


if __name__ == "__main__":
    main()

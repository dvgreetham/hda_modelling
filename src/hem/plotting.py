import pandas as pd
import numpy as np
import altair as alt


def plot_simulated_cost_results(costs_by_age, cost_names):
    cost_vars = ["total_cost_per_month"] + cost_names

    plot_data = (
        costs_by_age[["age"] + cost_vars]
        .groupby("age")
        .agg("mean")
        .reset_index()
        .melt(
            id_vars="age",
            value_vars=cost_vars,
            var_name="cost_item",
            value_name="cost_per_month",
        )
    )

    chart = (
        alt.Chart(plot_data)
        .mark_line()
        .encode(
            x="age",
            y="cost_per_month",
            color=alt.Color("cost_item", scale=alt.Scale(scheme="category20")),
        )
        .properties(width=500, height=400)
        .configure_axis(labelFontSize=20, titleFontSize=20)
        .configure_legend(
            labelFontSize=20,
            titleFontSize=20,
            labelLimit=0,
            symbolSize=300,
            symbolStrokeWidth=5,
        )
    )

    return chart, plot_data

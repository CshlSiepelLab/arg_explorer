def plot_horizontal_blocks(ax, df, label, color, metric, interval_start_col, interval_end_col):
    #sort on interval start
    df.sort_values(interval_start_col, inplace = True) 

    #plot horizontal lines
    for i, row in df.iterrows():
        ax.hlines(
            y=row[metric],
            xmin=row[interval_start_col],
            xmax=row[interval_end_col],
            color=color, lw=2,
            label=label if i == df.index[0] else None
        )

    # Vertical connectors
    for i in range(len(df) - 1):
        x_next = df.iloc[i + 1][interval_start_col]
        y1, y2 = df.iloc[i][metric], df.iloc[i + 1][metric]
        ax.vlines(x_next, ymin=min(y1, y2), ymax=max(y1, y2), color=color, lw=1.5)

##Specify exposed functions:
__all__ = [
    "plot_horizontal_blocks"
]
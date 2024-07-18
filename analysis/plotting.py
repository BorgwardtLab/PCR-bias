def standardize_plot(fig):
    fig.update_layout(
        template="simple_white",
        font_family="Inter",
        legend_font_size=28/3,
    )
    fig.update_yaxes(
        minor_ticks="outside", 
        title_font_family="Inter", 
        title_font_size=28/3, 
        tickfont_size=28/3, 
    )
    fig.update_xaxes(
        minor_ticks="outside", 
        title_font_family="Inter", 
        title_font_size=28/3, 
        tickfont_size=28/3, 
    )
    fig.for_each_annotation(lambda a: a.update(
        font_size=28/3,
        font_family="Inter",
    ))
    return fig
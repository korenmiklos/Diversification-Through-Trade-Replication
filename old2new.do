clear all
import delimited "output/sterr.csv"

foreach X in trade_barriers diversification {
    // Calculate upper and lower confidence intervals
    gen upper_ci = `X'_mean + 1.96 * `X'_se
    gen lower_ci = `X'_mean - 1.96 * `X'_se

    scatter `X'_mean `X', title(Percentage change in volatility due to `X')

    // Generate the graph
    twoway (scatter `X'_mean `X', mcolor(blue) msymbol(o)) ///
        (rcap upper_ci lower_ci `X', color(blue)) ///
        (line `X' `X', sort color(red)), ///
        legend(off) ///
        xtitle("Published version") ///
        ytitle("Julia 1.10 version") ///
        title("Percentage change in volatility due to `X'") ///
        subtitle("Error bars represent ±1.96 SE")

    // Create the graph
    graph export "output/`X'.png", width(1040) replace
    drop upper_ci lower_ci
}


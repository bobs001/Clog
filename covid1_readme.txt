################################################################
# How to calculate the mortality rate during an outbreak
#
# At present, it is tempting to estimate the case fatality rate by
# dividing the number of known deaths by the number of confirmed
# cases. The resulting number, however, does not represent the true
# case fatality rate and might be off by orders of magnitude [...]
#
# A precise estimate of the case fatality rate is therefore impossible
# at present.
#
# 2019-Novel Coronavirus (2019-nCoV): estimating the case
# fatality rate – a word of caution - Battegay Manue et al.,
# Swiss Med Wkly, February 7, 2020
#
# The case fatality rate (CFR) represents the proportion of cases who
# eventually die from a disease.
#
# Once an epidemic has ended, it is calculated with the formula:
# deaths/cases.
#
# But while an epidemic is still ongoing, as it is the case with the
# current novel coronavirus outbreak, this formula is, at the very
# least, "naïve" and can be "misleading if, at the time of analysis,
# the outcome is unknown for a non negligible proportion of patients."
# [8] (Methods for Estimating the Case Fatality Ratio for a Novel,
# Emerging Infectious Disease - Ghani et al, American Journal of
# Epidemiology).
#
# In other words, current deaths belong to a total case figure of the
# past, not to the current case figure in which the outcome (recovery
# or death) of a proportion (the most recent cases) hasn't yet been
# determined.
#
# The correct formula, therefore, would appear to be:
# 
# CFR = deaths at day.x / cases at day.x-{T}
# (where T = average time period from case confirmation to death)
#
# This would constitute a fair attempt to use values for cases
# and deaths belonging to the same group of patients.
#
# One issue can be that of determining whether there is enough data to
# estimate T with any precision, but it is certainly not T = 0 (what is
# implicitly used when applying the formula current deaths / current
# cases to determine CFR during an ongoing outbreak).
#
# Let's take, for example, the data at the end of February 8, 2020:
# 813 deaths (cumulative total) and 37,552 cases (cumulative total)
# worldwide.
#
# If we use the formula (deaths / cases) we get:
#
# 813 / 37,552 = 2.2% CFR (flawed formula).
#
# With a conservative estimate of T = 7 days as the average period from
# case confirmation to death, we would correct the above formula by
# using February 1 cumulative cases, which were 14,381, in the
# denominator:
#
# Feb. 8 deaths / Feb. 1 cases = 813 / 14,381 = 5.7% CFR (correct
# formula, and estimating T=7).
#
# T could be estimated by simply looking at the value of (current
# total deaths + current total recovered) and pair it with a case
# total in the past that has the same value. For the above formula,
# the matching dates would be January 26/27, providing an estimate for
# T of 12 to 13 days. This method of estimating T uses the same logic
# of the following method, and therefore will yield the same result.
#
# An alternative method, which has the advantage of not having to
# estimate a variable, and that is mentioned in the American Journal
# of Epidemiology study cited previously as a simple method that
# nevertheless could work reasonably well if the hazards of death and
# recovery at any time t measured from admission to the hospital,
# conditional on an event occurring at time t, are proportional, would
# be to use the formula:
#
# CFR = deaths / (deaths + recovered)
#
# which, with the latest data available, would be equal to:
#
# 15,408 / (15,408 + 100,605) = 13% CFR (worldwide)
#
# If we now exclude cases in mainland China, using current data on
# deaths and recovered cases, we get:
#
# 12,138 / (12,138 + 27,902) = 30.3% CFR (outside of mainland China)
#
# The sample size above is limited, and the data could be inaccurate
# (for example, the number of recoveries in countries outside of China
# could be lagging in our collection of data from numerous sources,
# whereas the number of cases and deaths is more readily available and
# therefore generally more up to par).
#
# There was a discrepancy in mortality rates (with a much higher
# mortality rate in China) which however is not being confirmed as the
# sample of cases outside of China is growing in size. On the
# contrary, it is now higher outside of China than within.
#
# That initial discrepancy was generally explained with a higher case
# detection rate outside of China especially with respect to Wuhan,
# where priority had to be initially placed on severe and critical
# cases, given the ongoing emergency.
#
# Unreported cases would have the effect of decreasing the denominator
# and inflating the CFR above its real value. For example, assuming
# 10,000 total unreported cases in Wuhan and adding them back to the
# formula, we would get a CFR of 12.2% (quite different from the CFR
# of 13% based strictly on confirmed cases).
#
# Neil Ferguson, a public health expert at Imperial College in the UK,
# said his “best guess” was that there were 100,000 affected by the
# virus even though there were only 2,000 confirmed cases at the
# time. [11]
#
# Without going that far, the possibility of a non negligible number
# of unreported cases in the initial stages of the crisis should be
# taken into account when trying to calculate the case fatally rate.
#
# As the days go by and the city organized its efforts and built the
# infrastructure, the ability to detect and confirm cases improved. As
# of February 3, for example, the novel coronavirus nucleic acid
# testing capability of Wuhan had increased to 4,196 samples per day
# from an initial 200 samples.[10]
#
# A significant discrepancy in case mortality rate can also be
# observed when comparing mortality rates as calculated and reported
# by China NHC: a CFR of 3.1% in the Hubei province (where Wuhan, with
# the vast majority of deaths is situated), and a CFR of 0.16% in
# other provinces (19 times less).
#
# Finally, we shall remember that while the 2003 SARS epidemic was
# still ongoing, the World Health Organization (WHO) reported a
# fatality rate of 4% (or as low as 3%), whereas the final case
# fatality rate ended up being 9.6%.

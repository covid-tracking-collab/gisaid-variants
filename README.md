# gisaid-variants

Repo containing work that merges basic GISAID metadata processing with Our World In Data cases
and population counts.

## Filtering GISAID sequences

These are the rules we're using to filter GISAID sequences:
- at least 20000 in length
- collection date must not be in the future and must be at the granularity of year/month/day
- earliest Dec 2019
- excluding sequences from the Nextstrain exclude list
- only human samples

Filters for future consideration (not doing these yet):
- excluding sequences with greater than 5% ambiguous base calls (N)

## Processed data column names and descriptions

Column name | In daily and/or weekly CSV | Description
------------|------| ------
`gisaid_collect_date` | daily | Same as "Collection date" value from original metadata file
`gisaid_country` | daily & weekly | 2nd string segment delimited by "/" in the "Location" value and manually changed to exactly match the `owid_location` value if there is a corresponding OWID location. If there is no match to OWID, this exists as a separate row for GISAID sequence metadata alone.
`All lineages` | daily & weekly | Count of all sequences regardless of lineage
`AY.1`, `AY.2`, ... `P.2`, `P.3` | daily & weekly | Variable number of break-out columns alphabetically ordered for sequence counts of each Pango lineage designated in the script as a VOC, VOI, or other important lineage
`Other lineages` | daily & weekly | `All lineages` minus the sum of counts across the `AY.1` ... `P.3` columns
`owid_date` | daily | Date of update from Our World In Data (OWID) file
`owid_new_cases` | daily & weekly | OWID's "new_cases". In daily file, this is the daily new cases for each `owid_date`. In weekly file, this is the week-by-week sum of new cases starting every Monday corresponding to `gisaid_collect_weekstartdate`.
`owid_new_cases_smoothed` | daily | "new_cases_smoothed" from OWID
`owid_location` | daily & weekly | OWID's "location". Should always match `gisaid_country`.  
`owid_continent` | daily & weekly | OWID's "continent".
`owid_population` | daily & weekly | OWID's "population".
`gisaid_collect_yearweek` | daily & weekly | the year and numbered week using ISO 8601 year and week (Monday as the first day of the week. Week 01 is the week containing Jan 4). I.e 2021-W01
`gisaid_collect_weekstartdate` | daily & weekly | the Monday date of each GISAID specimen collection week
`who_region` | daily & weekly | the designated WHO region for each `owid_location` according to https://ourworldindata.org/world-region-map-definitions
`aggregate_location` | daily & weekly | used for aggregate rows for continent- and global-level, blank for country-level rows
`gisaid_lagdays_min` | daily & weekly | precalculated summary stats for the length in days between specimen "Collection date" and sequence "Submission date" for all sequences per location binned daily or week-to-week. 
`gisaid_lagdays_q1` | daily & weekly | ^
`gisaid_lagdays_median` | daily & weekly | ^
`gisaid_lagdays_q3` | daily & weekly | ^
`gisaid_lagdays_max` | daily & weekly | ^
`owid_new_people_vaccinated` | daily & weekly | In daily file, this is calculated by taking the daily diff of "people_vaccinated" from OWID data. In weekly file, this is the week-by-week sum of the daily `owid_new_people_vaccinated` value per location.
`owid_new_people_fully_vaccinated` | daily & weekly | In daily file, this is calculated by taking the daily diff of "people_fully_vaccinated" from OWID data. In weekly file, this is the week-by-week sum of the daily `owid_new_people_fully_vaccinated` value per location.
`owid_people_fully_vaccinated` | weekly | The cumulative number of people fully vaccinated per location at the end of each week. This corresponds to the Sunday value reported by OWID for each country or continent.
`owid_people_fully_vaccinated_per_hundred` | weekly | `owid_people_fully_vaccinated` per 100 `owid_population` rounded to 2 decimal places.
`owid_people_vaccinated` | weekly | The cumulative number of people receiving at least one dose of vaccine per location at the end of each week. This corresponds to the Sunday value reported by OWID for each country or continent.
`owid_people_vaccinated_per_hundred`| weekly | `owid_people_vaccinated` per 100 `owid_population` rounded to 2 decimal places.
`sequences_over_new_cases` | weekly | `All lineages` divided by `owid_new_cases`
`new_cases_per_mil` | weekly | `owid_new_cases` per 1,000,000 `owid_population`
`who_alpha` | weekly | Sequence count for all Pango lineages designated by WHO under each greek-letter variant category at https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/ 
`who_beta` | weekly | ^
`who_gamma` | weekly | ^
`who_delta` | weekly | ^
`who_allvois` | weekly | Sequence count for all VOIs designated by WHO
`who_other` | weekly | `All lineages` minus all greek-lettered columns minus `who_allvois`

## Processing workflow overview

1. Load, clean up locations, and filter raw GISAID metadata file
    - load metadata.tsv from `--gisaid-metadata-file` filepath
    - optionally save this filtered sequence-level file locally with flag `--save-filtered-metadata`
2. Aggregate sequence counts to daily and country-level
    - break out Pango lineage columns designated as VOCs, VOIs, or other important lineages
3. Load OWID data and join to GISAID aggregate data based on country name and collection date
4. Aggregate to WHO region, continent, and global daily data and add as new `aggregate_location` rows
5. Calculate and add `gisaid_lagdays` daily summary stats by location
6. Save daily CSV to `--merged-gisaid-owid-out` filepath
7. If `--make-weekly-file` flag is set, continue to make weekly file
8. Aggregate all daily sequence, case, and vax counts to weekly counts based on `gisaid_collect_weekstartdate`
9. Overwrite continent- and global-level vax data with the toplines reported by OWID
10. Calculate and add `gisaid_lagdays` weekly summary stats by location
11. Calculate and add `sequences_over_new_cases` and `new_cases_per_mil`
12. Add `who_` greek-lettered variant columns
13. Saves to same path as `--merged-gisaid-owid-out` with "_weekly" appended to CSV filename.


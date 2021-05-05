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

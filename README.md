# Introducing Ecoshield Analytics: An App for Quantifying Control of Protected Areas
Natural resource management and conservation require effective monitoring tools to ensure the protection and sustainability of biodiversity within protected areas. We introduce Ecoshield Analytics, an innovative application designed to quantify the control of protected areas based on four types of indicators: patrol coverage, negative events, positive events, and spatially implicit indicators. These indicators are aggregated into an overall control score, providing a reproducible, data-driven monitoring tool that captures changes holistically. EcoShield Analytics addresses the needs of donors tracking conservation projects and the emerging biodiversity credit market requiring comprehensive state-of-protection indicators. This paper details the application's design, implementation, and potential impact on conservation efforts.

![image](https://github.com/Manuel-Weber-ETH/ecoshieldanalytics/assets/118481837/994029fa-208e-4d3d-8aef-9de22e277382)

This repository contains the script to run the app: EcoShieldAnalytics.py

It also contains example files to run the app:
- plateau.shp (and associated files) is a shapefile of a protected area
- positive_events.csv contains a series of coordinates of punctual events (type 2 and 3 indicators)
- patroltracks.kml contains a gps patrol track (type 1 indicator)

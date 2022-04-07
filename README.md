# Progenitor Search Tool

## Query input format 
Also see `sample_query.txt`

```
1               # 0: Neutron Star or 1: Black Hole
0.0,8.0         # Donor Mass Range (Msol)
6.0,14.0        # Accretor Mass Range (Msol)
-15.0,-4.0      # log10(MT Rate) Range (Msol/yr)
0.0,100000.0    # Orbital Period Range (days)
3.0,5.0         # log10(Donor Effective Temperature) Range (K)
```

Enter the quantities in order, if you don't have error estimates for a quantity, enter limits that span the entire range of possible values.

## Run from command line
Export the path to the track files as `DB_LOCATION` environment variable.
The path to query file is the first command line argument, as follows:

```console
 DB_LOCATION='/path/to/database/' ./progenlib.py ./sample_query.txt
```


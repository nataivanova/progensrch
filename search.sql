SELECT *
FROM tracks
WHERE mdonor    >= m1_lower   AND mdonor    <= m1_upper
AND   maccretor >= m2_lower   AND maccretor <= m2_upper
AND   mtrate    >= mt_lower   AND mtrate    <= mt_upper
AND   period    >= p_lower    AND period    <= p_upper
AND   teff      >= teff_lower AND teff      <= teff_upper;

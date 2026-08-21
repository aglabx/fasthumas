use crate::bed::Rgb;
use once_cell::sync::Lazy;
use regex::Regex;

/// A compiled color rule: regex pattern + RGB color.
struct ColorRule {
    re: Regex,
    color: Rgb,
}

// ---------------------------------------------------------------------------
// PROVENANCE: the color table below is NOT original to this repository. It is
// transplanted from hmmertblout2bed.awk by Lev I. Uralsky (Institute of
// Molecular Genetics, Moscow) -- https://github.com/enigene/hmmertblout2bed --
// including the later additions made in fedorrik/HumAS-HMMER_for_AnVIL.
// Neither upstream repository carries a license file, so the MIT grant in
// LICENSE does not extend to this table. See NOTICE.md before reusing it.
// ---------------------------------------------------------------------------

/// All 218+ color patterns from hmmertblout2bed.awk, in source order.
/// AWK uses IGNORECASE=1, so all patterns are compiled case-insensitive.
/// First match wins.
static COLOR_MAP: Lazy<Vec<ColorRule>> = Lazy::new(|| {
    let entries: Vec<(&str, u8, u8, u8)> = vec![
        ("M1", 255, 255, 0),
        ("Qa", 165, 39, 0),
        ("Pa", 75, 59, 162),
        ("Ta", 75, 59, 162),
        ("Ea", 153, 153, 255),
        ("Fa", 0, 128, 128),
        ("Ia", 153, 51, 102),
        ("Aa", 172, 172, 172),
        ("Ja", 225, 126, 231),
        ("Ba", 255, 153, 0),
        ("Ca", 224, 0, 64),
        ("Oa", 255, 229, 153),
        ("Na", 32, 160, 64),
        ("Ka", 0, 255, 0),
        ("Ha", 127, 96, 0),
        ("Ga", 255, 255, 0),
        ("R1", 0, 96, 192),
        ("R2", 102, 153, 255),
        ("D1", 153, 0, 255),
        ("D2", 210, 110, 250),
        ("D3", 255, 187, 255),
        ("D4", 205, 150, 205),
        ("D5", 181, 102, 195),
        ("D6", 159, 84, 183),
        ("D7", 176, 63, 176),
        ("D8", 141, 18, 141),
        ("D9", 102, 30, 102),
        ("FD", 139, 102, 139),
        ("W", 0, 255, 255),
        ("J1", 234, 153, 153),
        ("J2", 255, 204, 204),
        ("J3", 255, 20, 147),
        ("J4", 205, 16, 118),
        ("J5", 205, 96, 144),
        ("J6", 205, 140, 149),
        ("La", 128, 0, 128),
        (r"S4/6C13/14/21/22H2\.", 255, 10, 10),
        (r"S4/6C13/14/21H1\.", 255, 10, 10),
        (r"S4/6C13/14/21/22H8\.", 99, 171, 63),
        (r"S4C13/14/21/22H4\.", 94, 255, 148),
        (r"S4C13/14/21/22H5\.", 7, 171, 171),
        (r"S4C22H3\.", 7, 171, 171),
        (r"S4C13/14/21/22H9\.", 10, 92, 255),
        (r"S4C15H2\.", 120, 120, 171),
        (r"S4C15H2-B\.", 230, 179, 255),
        (r"S4C15H3\.", 230, 179, 255),
        (r"S4C15H3-A\.", 230, 179, 255),
        (r"S4C15H3-B\.", 171, 120, 137),
        (r"S4C15H3-AB\.", 180, 180, 180),
        (r"S4C20H4\.", 171, 63, 135),
        (r"S4C20H7\.", 255, 173, 10),
        (r"S4C20H8\.", 94, 255, 255),
        (r"S4C20H5-C\.", 7, 62, 171),
        (r"S4C20H5-D\.", 92, 10, 255),
        (r"S4CYH1L\.", 171, 7, 7),
        (r"S6C13/14/21/22H3\.", 255, 179, 230),
        (r"S6C13/14/21/22H3-A\.", 255, 179, 230),
        (r"S6C22H2-A\.", 255, 179, 230),
        (r"S6C13/14/21/22H3-B\.", 255, 94, 94),
        (r"S6C22H2-B\.", 255, 94, 94),
        (r"S5C11H3\.", 148, 94, 255),
        (r"S5C11H4\.", 148, 94, 255),
        (r"S5C5pH5\.", 171, 116, 7),
        (r"S5C5H5\.", 171, 116, 7),
        (r"S5C5/19pH5\.", 171, 116, 7),
        (r"S5C5/19H5\.", 171, 116, 7),
        (r"S5C5pH6\.", 173, 255, 10),
        (r"S5C5H6\.", 173, 255, 10),
        (r"S5C5pH7-A\.", 137, 171, 120),
        (r"S5C5H7-A\.", 137, 171, 120),
        (r"S5C5/19pH7-A\.", 137, 171, 120),
        (r"S5C5/19H7-A\.", 137, 171, 120),
        (r"S5C5pH7-B\.", 255, 179, 204),
        (r"S5C5H7-B\.", 255, 179, 204),
        (r"S5C5/19pH7-B\.", 255, 179, 204),
        (r"S5C5/19H7-B\.", 255, 179, 204),
        (r"S5C7H2\.", 178, 255, 204),
        (r"S5C5/19qH4-A\.", 63, 171, 171),
        (r"S5C5/19H4-A\.", 63, 171, 171),
        (r"S5C5/19qH4-B\.", 255, 10, 255),
        (r"S5C5/19H4-B\.", 255, 10, 255),
        (r"S5C13/14/21/22H6\.", 171, 63, 63),
        (r"S5C13/14/21H2\.", 171, 63, 63),
        (r"S5C20H6\.", 255, 201, 94),
        (r"S5C4H2\.", 135, 63, 171),
        (r"S5C5/19H9d\.", 255, 92, 10),
        (r"S5C1qH4d\.", 255, 10, 92),
        (r"S5C1H4d\.", 255, 10, 92),
        (r"S5C1qH5d\.", 10, 255, 10),
        (r"S5C1H5d\.", 10, 255, 10),
        (r"S5C1qH6d\.", 63, 99, 171),
        (r"S5C1H6d\.", 63, 99, 171),
        (r"S3C1pH2-B\.", 171, 7, 171),
        (r"S3C1H2-B\.", 171, 7, 171),
        (r"S3C1pH2-A\.", 62, 7, 171),
        (r"S3C1H2-A\.", 62, 7, 171),
        (r"S3C1qH2-D\.", 255, 179, 179),
        (r"S3C1H2-D\.", 255, 179, 179),
        (r"S3C1qH2-C\.", 171, 135, 63),
        (r"S3C1H2-C\.", 171, 135, 63),
        (r"S3C11H1L\.", 201, 255, 94),
        (r"S3CXH1L\.", 10, 255, 173),
        (r"S3C17H1-B\.", 120, 171, 171),
        (r"S3C17H1L\.", 179, 204, 255),
        (r"S3C17H1-C\.", 99, 63, 171),
        (r"S3C17H2\.", 255, 94, 255),
        (r"S3C17H2d\.", 255, 94, 255),
        (r"S3C11H3d\.", 171, 7, 62),
        (r"S2C2H1L\.", 171, 120, 120),
        (r"S2C2qH2-A\.", 255, 148, 94),
        (r"S2C2H2-A\.", 255, 148, 94),
        (r"S2C2pH2-B\.", 255, 230, 179),
        (r"S2C2H2-B\.", 255, 230, 179),
        (r"S2C4H1L\.", 135, 171, 63),
        (r"S2C8H1L\.", 94, 255, 94),
        (r"S2C9H1L\.", 7, 171, 116),
        (r"S2C15H1L\.", 10, 173, 255),
        (r"S2C13/21H1L\.", 120, 137, 171),
        (r"S2C13/21H1-B\.", 204, 179, 255),
        (r"S2C14/22H1L\.", 171, 63, 171),
        (r"S2C16pH2-B\.", 63, 171, 63),
        (r"S2C16H2-B\.", 63, 171, 63),
        (r"S2C16pH2-A\.", 255, 94, 148),
        (r"S2C16H2-A\.", 255, 94, 148),
        (r"S2C16pH3d\.", 171, 154, 120),
        (r"S2C16H3d\.", 171, 154, 120),
        (r"S2C16pH4d\.", 94, 255, 201),
        (r"S2C16H4d\.", 94, 255, 201),
        (r"S2CMH4d\.", 94, 255, 201),
        (r"S2C18pH2-A\.", 7, 116, 171),
        (r"S2C18H2-A\.", 7, 116, 171),
        (r"S2C18H1L\.", 10, 10, 255),
        (r"S2C18qH2-B\.", 255, 179, 255),
        (r"S2C18H2-B\.", 255, 179, 255),
        (r"S2C18qH2-D\.", 171, 63, 99),
        (r"S2C18H2-D\.", 171, 63, 99),
        (r"S2C18qH2-C\.", 94, 201, 255),
        (r"S2C18H2-C\.", 94, 201, 255),
        (r"S2C18qH2-E\.", 173, 10, 255),
        (r"S2C18H2-E\.", 173, 10, 255),
        (r"S2C20H2\.", 154, 171, 120),
        (r"S2C20H1L\.", 7, 7, 171),
        (r"S2C20H3\.", 178, 255, 178),
        (r"S02C20H3\.", 178, 255, 178),
        (r"S2CMH2d\.", 63, 171, 135),
        (r"S02CMH2d\.", 63, 171, 135),
        (r"S2C20H5d\.", 171, 120, 171),
        (r"S1C3H1L\.", 171, 171, 7),
        (r"S01/1C3H1L\.", 171, 171, 7),
        (r"S1C3H2\.", 92, 255, 10),
        (r"S01C3H2\.", 92, 255, 10),
        (r"S1C3H3d\.", 63, 63, 171),
        (r"S01C3/6H1d\.", 63, 63, 171),
        (r"S1CMH1d\.", 201, 94, 255),
        (r"S01CMH1d\.", 201, 94, 255),
        (r"S1C6H1L\.", 178, 255, 229),
        (r"S01C6H1L\.", 178, 255, 229),
        (r"S1C7H1L\.", 63, 135, 171),
        (r"S1C10H1L\.", 94, 94, 255),
        (r"S1C10H1-B\.", 116, 7, 171),
        (r"S1C10H1-C\.", 255, 10, 173),
        (r"S1C10H1-D\.", 171, 120, 137),
        (r"S1C10H2\.", 171, 99, 63),
        (r"S1C12H1L\.", 62, 171, 7),
        (r"S1C12H2\.", 10, 255, 92),
        (r"S1C16H1L\.", 255, 255, 10),
        (r"S1C1/5/19H1L\.", 120, 171, 154),
        (r"S1C5pH2\.", 179, 230, 255),
        (r"S1C5H2\.", 179, 230, 255),
        (r"S1C10/12H1d\.", 255, 204, 179),
        (r"S1CMH3d\.", 255, 204, 179),
        (r"S1C12H3d\.", 171, 7, 116),
        (r"S01C12H3d\.", 171, 7, 116),
        (r"S1C3H4-C\.", 10, 255, 255),
        // GREYR/YELLOW/YELLOWG/GREYB/GREYM/GREYZ variant patterns
        ("del_A_pos_70-GREYR.2", 212, 113, 0),
        ("ref_A_pos_70-GREYR.2", 14, 138, 30),
        ("del_27bp-GREYR.4", 87, 42, 212),
        ("ref_27bp-GREYR.4", 85, 212, 102),
        ("del_T_pos_150-YELLOW2.4", 106, 134, 212),
        ("ref_T_pos_150-YELLOW2.4", 83, 138, 119),
        ("del_7bp-YELLOW2.4", 148, 165, 212),
        ("ref_7bp-YELLOW2.4", 198, 212, 0),
        ("del_TTCT_pos_175-YELLOWG.4", 47, 14, 138),
        ("ref_TTCT_pos_175-YELLOWG.4", 212, 42, 155),
        ("ref_TC_pos_62-YELLOW2.4", 55, 138, 66),
        ("del_A_pos_63-GREYB.5", 69, 87, 138),
        ("ref_A_pos_63-GREYB.5", 127, 195, 212),
        ("alt_pos_10_69-GREYB.5", 165, 148, 212),
        ("ref_pos_10_69-GREYB.5", 85, 212, 0),
        ("del_C_pos_159-GREYM.5", 174, 21, 212),
        ("ref_C_pos_159-GREYM.5", 212, 42, 65),
        ("del_14bp_pos_49_62-GREYZ.5", 119, 85, 212),
        ("ref_14bp_pos_49_62-GREYZ.5", 87, 69, 138),
        ("del_T_pos_135-YELLOW2.6", 150, 127, 212),
        ("del_CCT_pos_60_62-YELLOW1.6", 107, 96, 138),
        ("ref_CCT_pos_60_62-YELLOW1.6", 0, 212, 141),
        ("del_AC_pos_60-GREYR.7", 113, 14, 138),
        ("ref_AC_pos_60-GREYR.7", 138, 93, 41),
        ("del_T_pos_135-YELLOW2.9", 121, 55, 138),
        ("ref_T_pos_135-YELLOW2.9", 134, 106, 212),
        ("del_A_pos_140-GREYP.10", 195, 127, 212),
        ("ref_A_pos_140-GREYP.10", 129, 96, 138),
        ("del_T_pos_135_and_pos_163-YELLOW2.10", 0, 110, 138),
        ("ref_T_pos_135_and_pos_163-YELLOW2.10", 138, 14, 96),
        ("del_T_pos_99-YELLOW3.11", 123, 212, 64),
        ("ref_T_pos_99-YELLOW3.11", 186, 85, 212),
        ("ins_SPL_W3J_GATC_pos_54-YELLOW2.10", 212, 106, 176),
        ("ref_SPL_W3J_GATC_pos_54-YELLOW2.10", 127, 83, 138),
        ("del_SPL_W4F_CCT_pos_61-YELLOW2.6", 199, 148, 212),
        ("ref_SPL_W4F_CCT_pos_61-YELLOW2.6", 0, 56, 212),
        ("alt_SPL_W4D_ATTTA_pos_59-YELLOW2.4", 138, 14, 30),
        ("ref_SPL_W4D_CTTTTCC_pos_59-YELLOW2.4", 64, 212, 162),
        // FEDOR updates
        (r"S3C1H2-F\.", 0, 255, 0),
        (r"S3C1H2-E\.", 0, 255, 255),
        // Note: S2С2H2-C has a Cyrillic С in the original AWK. We match both.
        (r"S2[CС]2H2-C\.", 0, 255, 0),
        (r"S3C11H2-A\.", 7, 171, 7),
        (r"S3C11H2-B\.", 0, 255, 255),
        (r"S4C15H4\.", 255, 150, 0),
    ];

    entries
        .into_iter()
        .map(|(pat, r, g, b)| {
            // Escape regex metacharacters in patterns that aren't already escaped.
            // The AWK patterns use literal `/` which isn't special in regex,
            // and `\.` which is already escaped.
            // We need case-insensitive matching.
            let re = regex::RegexBuilder::new(pat)
                .case_insensitive(true)
                .build()
                .unwrap_or_else(|e| panic!("bad color pattern '{}': {}", pat, e));
            ColorRule {
                re,
                color: Rgb(r, g, b),
            }
        })
        .collect()
});

/// Look up the color for a query name. First match wins.
/// Default: RGB(85,85,85) — grey.
pub fn get_color(query_name: &str) -> Rgb {
    for rule in COLOR_MAP.iter() {
        if rule.re.is_match(query_name) {
            return rule.color.clone();
        }
    }
    Rgb(85, 85, 85)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_known_colors() {
        assert_eq!(get_color("M1"), Rgb(255, 255, 0));
        assert_eq!(get_color("Qa"), Rgb(165, 39, 0));
        assert_eq!(get_color("S1C3H1L.1"), Rgb(171, 171, 7));
        assert_eq!(get_color("S4C15H2.3"), Rgb(120, 120, 171));
    }

    #[test]
    fn test_case_insensitive() {
        assert_eq!(get_color("m1"), Rgb(255, 255, 0));
        assert_eq!(get_color("qa"), Rgb(165, 39, 0));
    }

    #[test]
    fn test_default_color() {
        assert_eq!(get_color("ZZZZZ"), Rgb(85, 85, 85));
    }

    #[test]
    fn test_fedor_updates() {
        assert_eq!(get_color("S3C1H2-F.1"), Rgb(0, 255, 0));
        assert_eq!(get_color("S3C11H2-A.1"), Rgb(7, 171, 7));
        assert_eq!(get_color("S4C15H4.1"), Rgb(255, 150, 0));
    }
}

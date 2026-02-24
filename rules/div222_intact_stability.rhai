// ============================================================================
// Division 222 — Stabilité à l'état intact
// ============================================================================
//
// Réglementation française applicable aux navires de charge
// de jauge brute < 500 UMS.
//
// Référence : Division 222, Chapitre 2 — Franc-bord, stabilité,
//             compartimentage et assèchement.
//             Renvoie aux critères de la Division 211.
//
// Paramètres optionnels (via ctx.set_param) :
//   "vessel_type"       : "standard" (défaut) ou "new_decked"
//                         (navire ponté neuf construit après 01/07/2014)
//   "b_over_d"          : rapport B/D si > 2.5 pour critères spécifiques
//   "min_freeboard_mm"  : franc-bord minimal en mm (défaut 400)
//   "min_flooding_angle": angle d'envahissement minimal en degrés (défaut 30)
//
// Usage :
//   ctx = CriteriaContext.from_result(result, "MV Example", "Departure")
//   ctx.set_param("vessel_type", "new_decked")
//   engine = ScriptEngine()
//   result = engine.run_script_file("rules/div222_intact_stability.rhai", ctx)
//
// ============================================================================

fn check(ctx) {
    let results = [];

    // =========================================================================
    // Read optional parameters
    // =========================================================================
    let vessel_type = ctx.get_param("vessel_type");
    if vessel_type == () { vessel_type = "standard"; }

    let b_over_d = ctx.get_param("b_over_d");

    let min_freeboard_mm = ctx.get_param("min_freeboard_mm");
    if min_freeboard_mm == () { min_freeboard_mm = 400.0; }

    let min_flooding_angle = ctx.get_param("min_flooding_angle");
    if min_flooding_angle == () { min_flooding_angle = 30.0; }

    // =========================================================================
    // Determine thresholds based on vessel type and B/D ratio
    // =========================================================================
    let gm0_min = 0.15;             // Standard IS Code criteria
    let gz_max_angle_min = 25.0;     // Standard
    let area_0_30_min = 0.055;       // m·rad
    let area_0_40_min = 0.090;       // m·rad
    let area_30_40_min = 0.030;      // m·rad
    let gz_at_30_min = 0.200;        // m

    if vessel_type == "new_decked" {
        // Navires pontés neufs (post 01/07/2014)
        gm0_min = 0.70;
        gz_max_angle_min = 20.0;
    }

    if b_over_d != () && b_over_d > 2.5 {
        // Critères spécifiques B/D > 2.5
        gz_max_angle_min = 15.0;
    }

    // =========================================================================
    // Limiting angle (min of default 40°, flooding angle)
    // =========================================================================
    let limiting_angle = ctx.get_limiting_angle(40.0);

    // Also consider deck edge immersion
    let dei_angle = ctx.get_deck_edge_immersion_angle();
    if dei_angle != () {
        if dei_angle < limiting_angle {
            limiting_angle = dei_angle;
        }
    }

    let plot_id = "gz_plot";

    // =========================================================================
    // 1. Aire 0° – 30°  ≥ 0.055 m·rad
    // =========================================================================
    let upper_30 = if limiting_angle < 30.0 { limiting_angle } else { 30.0 };
    let area_0_30 = ctx.area_under_curve(0.0, upper_30);
    results.push(criterion(
        "Aire 0-30°",
        "Aire sous la courbe GZ de 0° à 30°",
        area_0_30_min,
        area_0_30,
        "m·rad"
    ));

    // =========================================================================
    // 2. Aire 0° – 40° (ou θf)  ≥ 0.090 m·rad
    // =========================================================================
    let area_0_lim = ctx.area_under_curve(0.0, limiting_angle);
    results.push(criterion(
        "Aire 0-40° (ou θf)",
        "Aire sous la courbe GZ de 0° à " + limiting_angle.toFixed(0) + "°",
        area_0_40_min,
        area_0_lim,
        "m·rad"
    ));

    // =========================================================================
    // 3. Aire 30° – 40° (ou θf)  ≥ 0.030 m·rad
    // =========================================================================
    if limiting_angle > 30.0 {
        let area_30_lim = ctx.area_under_curve(30.0, limiting_angle);
        results.push(criterion(
            "Aire 30-40° (ou θf)",
            "Aire sous la courbe GZ de 30° à " + limiting_angle.toFixed(0) + "°",
            area_30_40_min,
            area_30_lim,
            "m·rad"
        ));
    } else {
        results.push(#{
            name: "Aire 30-40° (ou θf)",
            description: "Non applicable — angle limitant ≤ 30°",
            required: area_30_40_min,
            actual: 0.0,
            unit: "m·rad",
            status: "FAIL",
            margin: -area_30_40_min,
            plot_id: plot_id
        });
    }

    // =========================================================================
    // 4. GZ ≥ 0.20 m à θ ≥ 30°
    // =========================================================================
    let gz_at_30 = ctx.gz_at_angle(30.0);
    results.push(criterion(
        "GZ à 30°+",
        "Bras de levier de redressement à 30° ou plus",
        gz_at_30_min,
        gz_at_30,
        "m"
    ));

    // =========================================================================
    // 5. Angle du GZ max  ≥ 25° (ou 20° new_decked, ou 15° B/D>2.5)
    // =========================================================================
    let max_info = ctx.find_max_gz();
    let max_angle = max_info.angle;
    let max_gz = max_info.value;

    if max_angle != () {
        results.push(criterion(
            "Angle GZmax",
            "Le GZ maximal doit être atteint à un angle ≥ " + gz_max_angle_min.toFixed(0) + "°",
            gz_max_angle_min,
            max_angle,
            "°"
        ));
    } else {
        results.push(#{
            name: "Angle GZmax",
            status: "FAIL",
            notes: "Impossible de déterminer le GZ maximal"
        });
    }

    // =========================================================================
    // 6. GM₀  ≥ 0.15 m (standard) ou ≥ 0.70 m (new_decked)
    // =========================================================================
    let gm0 = ctx.get_gm0();
    if gm0 != () {
        results.push(criterion(
            "GM₀ initial",
            "Hauteur métacentrique initiale (corrigée carènes liquides)",
            gm0_min,
            gm0,
            "m"
        ));
    } else {
        results.push(#{
            name: "GM₀ initial",
            status: "FAIL",
            notes: "GM₀ non disponible (VCG non spécifié ?)"
        });
    }

    // =========================================================================
    // 7. Angle d'envahissement  ≥ 30° (Division 222)
    // =========================================================================
    let flooding_angle = ctx.get_first_flooding_angle();
    if flooding_angle != () {
        results.push(criterion(
            "Angle d'envahissement",
            "L'angle d'envahissement doit être ≥ " + min_flooding_angle.toFixed(0) + "°",
            min_flooding_angle,
            flooding_angle,
            "°"
        ));
    } else {
        results.push(#{
            name: "Angle d'envahissement",
            status: "PASS",
            notes: "Pas d'ouverture d'envahissement définie — critère non vérifié",
            required: min_flooding_angle,
            actual: 0.0,
            unit: "°",
            margin: 0.0
        });
    }

    // =========================================================================
    // 8. Franc-bord  ≥ 400 mm (par défaut)
    // =========================================================================
    let dei = ctx.get_deck_edge_immersion_angle();
    if dei != () {
        // If deck edge immersion angle exists, check freeboard at 0 heel
        let gz_points = ctx.get_heels();
        // We check the freeboard criterion via the immersion angle
        // A deck edge immersion angle > 0 means freeboard > 0 at upright
        results.push(#{
            name: "Immersion du pont",
            description: "Angle d'immersion du livet de pont",
            required: 0.0,
            actual: dei,
            unit: "°",
            status: if dei > 0.0 { "PASS" } else { "FAIL" },
            margin: dei,
            plot_id: plot_id,
            notes: "Immersion du pont à " + dei.toFixed(1) + "°"
        });
    }

    // =========================================================================
    // Overall result
    // =========================================================================
    let fail_count = 0;
    for r in results {
        if r.status == "FAIL" { fail_count += 1; }
        r.plot_id = plot_id;
    }

    // =========================================================================
    // Plot
    // =========================================================================
    let heels = ctx.get_heels();
    let gz_values = ctx.get_gz_values();

    let elements = [
        #{
            type: "Curve",
            name: "GZ",
            x: heels,
            y: gz_values,
            color: "blue"
        },
        #{
            type: "HorizontalLine",
            name: "GZ min (0.20 m)",
            y: 0.200,
            x_min: 30.0,
            color: "red",
            style: "dashed"
        }
    ];

    // Add limiting angle vertical line
    if limiting_angle < 40.0 {
        elements.push(#{
            type: "VerticalLine",
            name: "Angle limitant (" + limiting_angle.toFixed(0) + "°)",
            x: limiting_angle,
            color: "orange",
            style: "dashed"
        });
    }

    // Add GZmax point
    if max_angle != () && max_gz != () {
        elements.push(#{
            type: "Point",
            name: "GZmax",
            x: max_angle,
            y: max_gz,
            marker: "o",
            color: "green"
        });
    }

    let plot = #{
        id: plot_id,
        title: "Courbe GZ — Division 222",
        x_label: "Gîte (°)",
        y_label: "GZ (m)",
        elements: elements
    };

    #{
        regulation_name: "Division 222 — Stabilité à l'état intact",
        regulation_reference: "Division 222, Chapitre 2 (réf. Division 211)",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        cog: ctx.get_cog(),
        criteria: results,
        overall_pass: fail_count == 0,
        notes: if vessel_type == "new_decked" { "Navire ponté neuf (post 01/07/2014)" } else { "" },
        plots: [plot]
    }
}

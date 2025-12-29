package main

import "core:fmt"
import "core:strings"
import "core:math"
import "core:time"
import "core:os"

DEG2RAD :: math.PI / 180.0

Vec2 :: struct {

    x : f64,
    y : f64
}

Planet :: struct {

	name      : string,
	color     : string,
	radius_px : f64,

	// Keplerian elements at J2000 and rates per century
	a0        : f64,
    a1        : f64,
	e0        : f64,
    e1        : f64,
	I0        : f64,    // deg
    I1        : f64,    // deg
	L0        : f64,    // deg
    L1        : f64,    // deg
	peri0     : f64,    // longitude of perihelion, deg
    peri1     : f64,    // longitude of perihelion, deg
	node0     : f64,    // longitude of ascending node, deg
    node1     : f64,    // longitude of ascending node, deg

	// Extra mean-anomaly terms ( used for outer planets; inner planets = 0 )
	b         : f64,
    c         : f64,
    s         : f64,
    f         : f64,
}

wrap_deg :: proc ( x : f64 ) ->
                   f64 {

    y := math.mod( x, 360.0 )
	if y < 0 {

	    y += 360.0
	}

	return y
}

wrap_pm180 :: proc ( x : f64 ) ->
                     f64 {

	y := wrap_deg( x )
	if y > 180.0 {

	    y -= 360.0
	}

	return y
}

wrap_rad_pm_pi :: proc ( x : f64 ) ->
                         f64 {

	y := math.mod( x, 2.0 * math.PI )
	if y >  math.PI {

	    y -= 2.0 * math.PI
	}
	if y < -math.PI {

	    y += 2.0 * math.PI
	}

	return y
}

julian_day_ymd :: proc ( year  : int,
                         month : int,
                         day   : int ) ->
                         f64 {

	// Gregorian calendar, 0h UT
	y := year
	m := month
	d := f64( day )

	if m <= 2 {

		y -= 1
		m += 12
	}

	A := y / 100
	B := 2 - A + ( A / 4 )

	jd1 := f64( i64( 365.25 * f64( y + 4716 ) ) )
	jd2 := f64( i64( 30.6001 * f64( m + 1 ) ) )

	return jd1 + jd2 + d + f64( B ) - 1524.5
}

solve_kepler_E :: proc ( M : f64,
                         e : f64 ) ->
                         f64 {

	E := M
	for _ in 0 ..< 10 {

		f  := E - e * math.sin( E ) - M
		fp := 1.0 - e * math.cos( E )
		E  -= f / fp
	}

	return E
}

// Radial log mapping : ( x, y ) in AU -> ( x, y ) in px
// r' = ln( 1 + r ) / ln( 1 + r_max ) * r_px_max, angle preserved
map_log_radial :: proc( x_au     : f64,
                        y_au     : f64,
                        r_max_au : f64,
                        r_px_max : f64 ) ->
                        Vec2 {

	r := math.sqrt( x_au * x_au + y_au * y_au )
	if r <= 1e-12 {

		return Vec2{ 0, 0 }
	}
	rn := math.ln( 1.0 + r ) / math.ln( 1.0 + r_max_au )      // 0 .. 1
	rp := rn * r_px_max
	return Vec2{ ( x_au / r ) * rp, ( y_au / r ) * rp }
}

// Heliocentric ecliptic coordinates in AU ( x, y ) using Keplerian elements
planet_xy_au :: proc ( p  : Planet,
                       jd : f64 ) ->
                       Vec2 {

	T := ( jd - 2451545.0 ) / 36525.0

	a    := p.a0    + p.a1 * T
	e    := p.e0    + p.e1 * T
	Ideg := p.I0    + p.I1 * T
	Ldeg := p.L0    + p.L1 * T
	peri := p.peri0 + p.peri1 * T
	node := p.node0 + p.node1 * T

	Ideg = wrap_deg( Ideg )
	Ldeg = wrap_deg( Ldeg )
	peri = wrap_deg( peri )
	node = wrap_deg( node )

	// Mean anomaly ( deg ), with optional outer-planet corrections
	Mdeg := ( Ldeg - peri )
	if p.b != 0 || p.c != 0 || p.s != 0 || p.f != 0 {

		Mdeg += p.b * T * T + p.c * math.cos( ( p.f * T ) * DEG2RAD ) + p.s * math.sin( ( p.f * T ) * DEG2RAD )
	}

	Mdeg = wrap_pm180( Mdeg )
	M := wrap_rad_pm_pi( Mdeg * DEG2RAD )

	E := solve_kepler_E( M, e )

	xp := a * ( math.cos( E ) - e )
	yp := a * math.sqrt( 1.0 - e * e ) * math.sin( E )

	omega := wrap_deg( peri - node ) * DEG2RAD
	Omega := node * DEG2RAD
	I     := Ideg * DEG2RAD

	cw := math.cos( omega )
	sw := math.sin( omega )
	cO := math.cos( Omega )
	sO := math.sin( Omega )
	cI := math.cos( I )

	x := ( cw * cO - sw * sO * cI ) * xp + ( -sw * cO - cw * sO * cI ) * yp
	y := ( cw * sO + sw * cO * cI ) * xp + ( -sw * sO + cw * cO * cI ) * yp

	return Vec2{ x, y }
}

compute_rmax_au :: proc ( planets : [ ]Planet,
                          Tmid    : f64 ) ->
                          f64 {

	rmax := 0.0
	for p in planets {
		a := p.a0 + p.a1 * Tmid
		e := p.e0 + p.e1 * Tmid
		rap := a * ( 1.0 + e )     // aphelion distance ( AU )
		if rap > rmax {

		    rmax = rap
		}
	}

	// Small padding so outer orbit doesn't kiss the edge
	return rmax * 1.03
}

write_orbit_path_log :: proc ( b        : ^strings.Builder,
                               p        : Planet,
                               Tmid     : f64,
                               r_max_au : f64,
                               r_px_max : f64 ) {

	a    := p.a0    + p.a1 * Tmid
	e    := p.e0    + p.e1 * Tmid
	Ideg := wrap_deg( p.I0    + p.I1 * Tmid )
	peri := wrap_deg( p.peri0 + p.peri1 * Tmid )
	node := wrap_deg( p.node0 + p.node1 * Tmid )

	omega := wrap_deg( peri - node ) * DEG2RAD
	Omega := node * DEG2RAD
	I     := Ideg * DEG2RAD

	cw := math.cos( omega )
	sw := math.sin( omega )
	cO := math.cos( Omega )
	sO := math.sin( Omega )
	cI := math.cos( I )

	steps := 720
	fmt.sbprintf( b, "<path fill=\"none\" stroke=\"%s\" opacity=\"0.42\" stroke-width=\"1.6\" vector-effect=\"non-scaling-stroke\" d=\"",
		          p.color )

	for i in 0 ..= steps {

		v := ( 2.0 * math.PI ) * ( f64( i ) / f64( steps ) )
		r := a * ( 1.0 - e * e ) / ( 1.0 + e * math.cos( v ) )

		xp := r * math.cos( v )
		yp := r * math.sin( v )

		x := ( cw * cO - sw * sO * cI ) * xp + ( -sw * cO - cw * sO * cI ) * yp
		y := ( cw * sO + sw * cO * cI ) * xp + ( -sw * sO + cw * cO * cI ) * yp

		m  := map_log_radial( x, y, r_max_au, r_px_max )
		sx := m.x
		sy := -m.y // Flip y for screen coords

		if i == 0 {

			fmt.sbprintf( b, "M %.1f %.1f", sx, sy )
		} else {

			fmt.sbprintf( b, " L %.1f %.1f", sx, sy )
		}
	}

	strings.write_string( b, "\"/>\n" )
}

write_animate_values_log :: proc (
                           	b         : ^strings.Builder,
                           	attr      : string,
                           	p         : Planet,
                           	jd_start  : f64,
                            jd_end    : f64,
                           	step_days : f64,
                           	r_max_au  : f64,
                            r_px_max  : f64,
                           	dur_s     : int ) {

	span := jd_end - jd_start
	n := int( math.floor( span / step_days ) ) + 1
	if n < 2 {

	    n = 2
	}

	fmt.sbprintf( b, "<animate attributeName=\"%s\" dur=\"%ds\" repeatCount=\"indefinite\" calcMode=\"linear\" values=\"",
		          attr, dur_s )

	for i in 0 ..< n {

		jd := jd_start + f64( i ) * step_days
		if jd > jd_end {

		    jd = jd_end
	    }

		pos := planet_xy_au( p, jd )
		m   := map_log_radial( pos.x, pos.y, r_max_au, r_px_max )

		val := m.x
		if attr == "cy" {

			val = -m.y
		}

		if i != 0 {

			strings.write_string( b, ";" )
		}

		fmt.sbprintf( b, "%.1f", val )
	}

	strings.write_string( b, "\"/>\n" )
}

main :: proc ( ) {

	// Layout
	VIEW_HALF  := 520.0
	R_PX_MAX   := 470.0    // max radius after log mapping
	STEP_DAYS  := 3.0      // smaller = smoother, bigger SVG
	DUR_S      := 90

	// Time range: 1976-01-01 -> today
	now      := time.now( )
	y, m, d  := time.date( now )
	jd_start := julian_day_ymd( 1976, 1, 1 )
	jd_end   := julian_day_ymd( y, int( m ), d )
	jd_mid   := 0.5 * ( jd_start + jd_end )
	Tmid     := ( jd_mid - 2451545.0 ) / 36525.0

	planets := [ ]Planet{

		{ "Mercury", "#b9c0c7", 3.2, 0.38709843,  0.00000000, 0.20563661,  0.00002123, 7.00559432, -0.00590158, 252.25166724, 149472.67486623,  77.45771895, 0.15940013,  48.33961819, -0.12214182, 0, 0, 0, 0 },
		{ "Venus",   "#ffd166", 4.0, 0.72332102, -0.00000026, 0.00676399, -0.00005107, 3.39777545,  0.00043494, 181.97970850,  58517.81560260, 131.76755713, 0.05679648,  76.67261496, -0.27274174, 0, 0, 0, 0 },
		{ "Earth",   "#4cc9f0", 4.2, 1.00000018, -0.00000003, 0.01673163, -0.00003661,-0.00054346, -0.01337178, 100.46691572,  35999.37306329, 102.93005885, 0.31795260,  -5.11260389, -0.24123856, 0, 0, 0, 0 },
		{ "Mars",    "#ff5d6c", 3.7, 1.52371243,  0.00000097, 0.09336511,  0.00009149, 1.85181869, -0.00724757,  -4.56813164,  19140.29934243, -23.91744784, 0.45223625,  49.71320984, -0.26852431, 0, 0, 0, 0 },
		{ "Jupiter", "#f4a261", 6.2, 5.20248019, -0.00002864, 0.04853590,  0.00018026, 1.29861416, -0.00322699,  34.33479152,   3034.90371757,  14.27495244, 0.18199196, 100.29282654,  0.13024619, -0.00012452, 0.06064060, -0.35635438, 38.35125000 },
		{ "Saturn",  "#e9c46a", 5.6, 9.54149883, -0.00003065, 0.05550825, -0.00032044, 2.49424102,  0.00451969,  50.07571329,   1222.11494724,  92.86136063, 0.54179478, 113.63998702, -0.25015002,  0.00025899,-0.13434469,  0.87320147, 38.35125000 },
		{ "Uranus",  "#80ffdb", 5.0,19.18797948, -0.00020455, 0.04685740, -0.00001550, 0.77298127, -0.00180155, 314.20276625,    428.49512595, 172.43404441, 0.09266985,  73.96250215,  0.05739699,  0.00058331,-0.97731848,  0.17689245,  7.67025000 },
		{ "Neptune", "#7c4dff", 5.0,30.06952752,  0.00006447, 0.00895439,  0.00000818, 1.77005520,  0.00022400, 304.22289287,    218.46515314,  46.68158724, 0.01009938, 131.78635853, -0.00606302, -0.00041348, 0.68346318, -0.10162547,  7.67025000 },
	}

	r_max_au := compute_rmax_au( planets, Tmid )

	b := strings.builder_make( )
	defer strings.builder_destroy( & b )

	fmt.sbprintf( & b, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"-%.0f -%.0f %.0f %.0f\" width=\"500\" height=\"500\">\n",
	              VIEW_HALF, VIEW_HALF, 2 * VIEW_HALF, 2 * VIEW_HALF )

	strings.write_string( & b,
		"<defs>\n" +
		"  <filter id=\"glow\" x=\"-40%\" y=\"-40%\" width=\"180%\" height=\"180%\">\n" +
		"    <feGaussianBlur stdDeviation=\"2\" result=\"b\"/>\n" +
		"    <feMerge><feMergeNode in=\"b\"/><feMergeNode in=\"SourceGraphic\"/></feMerge>\n" +
		"  </filter>\n" +
		"</defs>\n" )

	// Black background
	strings.write_string( & b, "<rect x=\"-520\" y=\"-520\" width=\"1040\" height=\"1040\" rx=\"24\" fill=\"#000000\"/>\n" )

	// Orbits (log-scaled)
	strings.write_string( & b, "<g id=\"orbits\">\n" )
	for p in planets {

		write_orbit_path_log( & b, p, Tmid, r_max_au, R_PX_MAX )
	}
	strings.write_string( & b, "</g>\n" )

	// Sun
	strings.write_string( & b,
		"<g id=\"sun\">\n" +
		"  <circle cx=\"0\" cy=\"0\" r=\"12\" fill=\"#fff1a8\" filter=\"url(#glow)\"/>\n" +
		"  <circle cx=\"0\" cy=\"0\" r=\"30\" fill=\"none\" stroke=\"#fff1a8\" opacity=\"0.12\"/>\n" +
		"</g>\n" )

	// Planets ( animated, log-scaled )
	strings.write_string( & b, "<g id=\"planets\">\n" )
	for p in planets {

		start_pos := planet_xy_au( p, jd_start )
		mpos := map_log_radial( start_pos.x, start_pos.y, r_max_au, R_PX_MAX )
		cx0  := mpos.x
		cy0  := -mpos.y

		fmt.sbprintf( & b, "  <circle r=\"%.1f\" cx=\"%.1f\" cy=\"%.1f\" fill=\"%s\" filter=\"url(#glow)\">\n",
			          p.radius_px, cx0, cy0, p.color )

		write_animate_values_log( & b, "cx", p, jd_start, jd_end, STEP_DAYS, r_max_au, R_PX_MAX, DUR_S )
		write_animate_values_log( & b, "cy", p, jd_start, jd_end, STEP_DAYS, r_max_au, R_PX_MAX, DUR_S )

		strings.write_string( & b, "  </circle>\n" )
	}
	strings.write_string( & b, "</g>\n" )

	strings.write_string( & b, "</svg>\n" )

	svg := strings.to_string( b )
	ok  := os.write_entire_file( "solar_orbits_log.svg", transmute( [ ]byte )( svg ) )
	if !ok {

		fmt.println( "Failed to write solar_orbits_log.svg" )
		return
	}
	fmt.println( "Wrote solar_orbits_log.svg" )
}

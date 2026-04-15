"""
caliper/viewer.py
Interactive 3D WebGL Inspection Dashboard Generator.

Generates a standalone, zero-dependency HTML5/WebGL inspection dashboard
for visual metrology analysis. Integrates Three.js for interactive 3D rendering,
vertex-colored deviation heatmaps, cross-sectional clipping planes, point probing,
and canvas-based deviation distribution histograms.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from caliper.pipeline import InspectionResult


def generate_html_viewer(
    result: InspectionResult,
    title: Optional[str] = None,
    tolerance_band: float = 0.5,
    max_display_points: int = 50000,
) -> str:
    """Generate self-contained HTML/WebGL inspection report.

    Args:
        result: The InspectionResult containing points, deviations, mesh, and GD&T outcomes.
        title: Optional custom dashboard title.
        tolerance_band: Bilateral tolerance band in millimeters for heatmap normalization.
        max_display_points: Point cloud subsampling limit for smooth browser rendering.

    Returns:
        Full HTML string containing the interactive 3D viewer.
    """
    pts = result.scanned_points.points
    devs = result.signed_deviations
    n_pts = len(pts)

    step = 1
    if n_pts > max_display_points:
        step = int(math.ceil(n_pts / float(max_display_points)))

    sampled_points = pts[::step]
    sampled_devs = devs[::step]

    point_coords_flat: List[float] = []
    for p in sampled_points:
        point_coords_flat.extend([float(p[0]), float(p[1]), float(p[2])])

    point_devs_list: List[float] = [float(d) for d in sampled_devs]

    mesh_vertices_flat: List[float] = []
    mesh_triangles_flat: List[int] = []
    if result.nominal_mesh is not None:
        for v in result.nominal_mesh.vertices:
            mesh_vertices_flat.extend([float(v[0]), float(v[1]), float(v[2])])
        for tri in result.nominal_mesh.triangle_indices:
            mesh_triangles_flat.extend([int(tri[0]), int(tri[1]), int(tri[2])])

    disposition_color = "#10b981" if result.disposition == "PASS" else ("#f59e0b" if result.disposition == "REWORK" else "#ef4444")
    dash_title = title or f"Metrology Inspection Dashboard - Part #{result.metadata.get('part_id', 'MET-001')}"

    gdt_rows_html = ""
    for gdt in result.gdt_results:
        status_badge = '<span class="badge pass">PASS</span>' if gdt.passed else '<span class="badge fail">FAIL</span>'
        gdt_rows_html += f"""
        <tr>
            <td><strong>{gdt.feature_name}</strong></td>
            <td>{gdt.characteristic}</td>
            <td>{gdt.measured_value:.4f} mm</td>
            <td>{gdt.tolerance_limit:.4f} mm</td>
            <td>{gdt.bonus_tolerance:.4f} mm</td>
            <td>{status_badge}</td>
        </tr>
        """

    defect_rows_html = ""
    if result.defects:
        for d in result.defects:
            defect_rows_html += f"""
            <tr>
                <td>#{d.cluster_id}</td>
                <td><span class="defect-type">{d.defect_type}</span></td>
                <td>{d.point_count}</td>
                <td>{d.peak_deviation:.4f} mm</td>
                <td>({d.bounding_box_min[0]:.1f}, {d.bounding_box_min[1]:.1f}, {d.bounding_box_min[2]:.1f})</td>
            </tr>
            """
    else:
        defect_rows_html = "<tr><td colspan='5' style='text-align:center; color:#9ca3af;'>No critical defect clusters detected</td></tr>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{dash_title}</title>
    <style>
        :root {{
            --bg-primary: #0a0e17;
            --bg-secondary: #111827;
            --bg-card: #1f2937;
            --text-main: #f3f4f6;
            --text-muted: #9ca3af;
            --border: #374151;
            --accent: #3b82f6;
            --pass: #10b981;
            --fail: #ef4444;
            --rework: #f59e0b;
        }}
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        }}
        body {{
            background-color: var(--bg-primary);
            color: var(--text-main);
            overflow: hidden;
            display: flex;
            height: 100vh;
            width: 100vw;
        }}
        #viewport-container {{
            flex: 1;
            position: relative;
            background: radial-gradient(circle at center, #1e293b 0%, #0f172a 100%);
            height: 100%;
        }}
        #canvas3d {{
            width: 100%;
            height: 100%;
            display: block;
        }}
        #sidebar {{
            width: 440px;
            background-color: var(--bg-secondary);
            border-left: 1px solid var(--border);
            overflow-y: auto;
            display: flex;
            flex-direction: column;
            gap: 16px;
            padding: 20px;
            z-index: 10;
        }}
        .header-card {{
            background: linear-gradient(135deg, #1e293b, #0f172a);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 16px;
        }}
        .header-title {{
            font-size: 1.15rem;
            font-weight: 700;
            margin-bottom: 6px;
            letter-spacing: -0.01em;
        }}
        .header-meta {{
            font-size: 0.82rem;
            color: var(--text-muted);
            display: flex;
            justify-content: space-between;
        }}
        .disposition-banner {{
            display: inline-block;
            margin-top: 10px;
            padding: 4px 12px;
            border-radius: 4px;
            font-weight: 700;
            font-size: 0.88rem;
            letter-spacing: 0.05em;
            background-color: {disposition_color}22;
            color: {disposition_color};
            border: 1px solid {disposition_color};
        }}
        .section-card {{
            background-color: var(--bg-card);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 14px;
        }}
        .section-title {{
            font-size: 0.88rem;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: var(--text-muted);
            margin-bottom: 12px;
            font-weight: 600;
            display: flex;
            justify-content: space-between;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 10px;
        }}
        .metric-box {{
            background: #111827;
            padding: 10px;
            border-radius: 6px;
            border: 1px solid #2d3748;
        }}
        .metric-label {{
            font-size: 0.72rem;
            color: var(--text-muted);
            text-transform: uppercase;
        }}
        .metric-value {{
            font-size: 1.1rem;
            font-weight: 700;
            margin-top: 2px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.8rem;
            margin-top: 4px;
        }}
        th, td {{
            text-align: left;
            padding: 6px 8px;
            border-bottom: 1px solid #2d3748;
        }}
        th {{
            color: var(--text-muted);
            font-weight: 600;
        }}
        .badge {{
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.72rem;
            font-weight: 700;
        }}
        .badge.pass {{ background-color: #065f46; color: #34d399; }}
        .badge.fail {{ background-color: #7f1d1d; color: #f87171; }}
        .defect-type {{
            background-color: #4b5563;
            padding: 2px 5px;
            border-radius: 3px;
            font-size: 0.72rem;
        }}
        .control-group {{
            margin-bottom: 10px;
        }}
        .control-label {{
            display: flex;
            justify-content: space-between;
            font-size: 0.78rem;
            color: var(--text-muted);
            margin-bottom: 4px;
        }}
        input[type="range"] {{
            width: 100%;
            height: 6px;
            background: #374151;
            border-radius: 3px;
            outline: none;
        }}
        .checkbox-label {{
            font-size: 0.8rem;
            display: flex;
            align-items: center;
            gap: 8px;
            cursor: pointer;
            margin-bottom: 6px;
        }}
        #color-legend {{
            position: absolute;
            bottom: 24px;
            left: 24px;
            background: rgba(17, 24, 39, 0.85);
            backdrop-filter: blur(8px);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 12px 16px;
            display: flex;
            flex-direction: column;
            gap: 6px;
            z-index: 5;
            width: 240px;
        }}
        .legend-gradient {{
            height: 12px;
            border-radius: 4px;
            background: linear-gradient(to right, #0000ff, #00ffff, #00ff00, #ffff00, #ff0000);
        }}
        .legend-labels {{
            display: flex;
            justify-content: space-between;
            font-size: 0.75rem;
            color: var(--text-muted);
        }}
        #probe-card {{
            position: absolute;
            top: 24px;
            left: 24px;
            background: rgba(17, 24, 39, 0.9);
            backdrop-filter: blur(8px);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 12px;
            display: none;
            z-index: 5;
            font-size: 0.8rem;
            min-width: 180px;
        }}
        #histogram-canvas {{
            width: 100%;
            height: 90px;
            background: #111827;
            border-radius: 4px;
            margin-top: 6px;
        }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="viewport-container">
        <canvas id="canvas3d"></canvas>
        <div id="color-legend">
            <span style="font-size: 0.75rem; font-weight: 600;">Signed Deviation Spectrum (mm)</span>
            <div class="legend-gradient"></div>
            <div class="legend-labels">
                <span id="legend-min">-{tolerance_band:.2f}</span>
                <span>0.00</span>
                <span id="legend-max">+{tolerance_band:.2f}</span>
            </div>
        </div>
        <div id="probe-card">
            <div style="font-weight:700; margin-bottom:4px; color:var(--accent);">Surface Point Probe</div>
            <div>X: <span id="probe-x">0.00</span> mm</div>
            <div>Y: <span id="probe-y">0.00</span> mm</div>
            <div>Z: <span id="probe-z">0.00</span> mm</div>
            <div>Dev: <span id="probe-dev" style="font-weight:700;">0.00</span> mm</div>
        </div>
    </div>

    <div id="sidebar">
        <div class="header-card">
            <div class="header-title">{dash_title}</div>
            <div class="header-meta">
                <span>ISO 1101 / ASME Y14.5</span>
                <span>{result.metadata.get("operator", "Automated Inspection")}</span>
            </div>
            <div class="disposition-banner">DISPOSITION: {result.disposition}</div>
        </div>

        <div class="section-card">
            <div class="section-title">Surface Error Statistics</div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">Points Evaluated</div>
                    <div class="metric-value">{result.statistics.point_count:,}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">RMSE Deviation</div>
                    <div class="metric-value">{result.statistics.rmse:.4f} mm</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Mean Deviation</div>
                    <div class="metric-value">{result.statistics.mean_deviation:+.4f} mm</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Std Deviation</div>
                    <div class="metric-value">{result.statistics.std_deviation:.4f} mm</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Max Positive (Burr)</div>
                    <div class="metric-value" style="color:#ef4444;">+{result.statistics.max_positive:.4f} mm</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Max Negative (Void)</div>
                    <div class="metric-value" style="color:#3b82f6;">{result.statistics.max_negative:.4f} mm</div>
                </div>
            </div>
            <canvas id="histogram-canvas"></canvas>
        </div>

        <div class="section-card">
            <div class="section-title">Display Controls</div>
            <label class="checkbox-label">
                <input type="checkbox" id="toggle-cad" checked> Show Nominal CAD Surface
            </label>
            <label class="checkbox-label">
                <input type="checkbox" id="toggle-wireframe"> CAD Wireframe Mode
            </label>
            <label class="checkbox-label">
                <input type="checkbox" id="toggle-points" checked> Show Scan Points
            </label>
            <div class="control-group" style="margin-top: 10px;">
                <div class="control-label">
                    <span>Point Size</span>
                    <span id="point-size-val">2.5 px</span>
                </div>
                <input type="range" id="slider-point-size" min="1" max="10" step="0.5" value="2.5">
            </div>
            <div class="control-group">
                <div class="control-label">
                    <span>Cross-Section Clip Z</span>
                    <span id="clip-z-val">Off</span>
                </div>
                <input type="range" id="slider-clip-z" min="-100" max="100" step="1" value="100">
            </div>
        </div>

        <div class="section-card">
            <div class="section-title">GD&T Tolerance Verification</div>
            <table>
                <thead>
                    <tr>
                        <th>Feature</th>
                        <th>Type</th>
                        <th>Meas</th>
                        <th>Tol</th>
                        <th>Bonus</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
                    {gdt_rows_html}
                </tbody>
            </table>
        </div>

        <div class="section-card">
            <div class="section-title">Morphological Defects ({len(result.defects)})</div>
            <table>
                <thead>
                    <tr>
                        <th>ID</th>
                        <th>Classification</th>
                        <th>Pts</th>
                        <th>Peak</th>
                        <th>BBox Origin</th>
                    </tr>
                </thead>
                <tbody>
                    {defect_rows_html}
                </tbody>
            </table>
        </div>
    </div>

    <script>
        const POINT_COORDS = {json.dumps(point_coords_flat)};
        const POINT_DEVS = {json.dumps(point_devs_list)};
        const MESH_VERTS = {json.dumps(mesh_vertices_flat)};
        const MESH_TRIS = {json.dumps(mesh_triangles_flat)};
        const TOLERANCE_BAND = {tolerance_band};

        const container = document.getElementById('viewport-container');
        const canvas = document.getElementById('canvas3d');

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(45, container.clientWidth / container.clientHeight, 0.1, 2000);
        camera.position.set(120, 100, 150);

        const renderer = new THREE.WebGLRenderer({{ canvas: canvas, antialias: true }});
        renderer.setSize(container.clientWidth, container.clientHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.localClippingEnabled = true;

        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;

        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        scene.add(ambientLight);
        const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight1.position.set(200, 300, 200);
        scene.add(dirLight1);
        const dirLight2 = new THREE.DirectionalLight(0x90b0d0, 0.4);
        dirLight2.position.set(-200, -100, -200);
        scene.add(dirLight2);

        const grid = new THREE.GridHelper(200, 40, 0x374151, 0x1f2937);
        grid.position.y = -40;
        scene.add(grid);

        const axes = new THREE.AxesHelper(30);
        axes.position.set(-80, -39, -80);
        scene.add(axes);

        const clipPlaneZ = new THREE.Plane(new THREE.Vector3(0, 0, -1), 1000);

        function turboColormap(t) {{
            t = Math.max(0, Math.min(1, t));
            let r = Math.sin(t * Math.PI - Math.PI / 2) * 0.5 + 0.5;
            let g = Math.sin(t * Math.PI) * 0.8;
            let b = Math.sin(t * Math.PI + Math.PI / 2) * 0.5 + 0.5;
            if (t < 0.25) {{
                return [0.1, 0.2 + 0.8 * (t / 0.25), 1.0];
            }} else if (t < 0.5) {{
                let u = (t - 0.25) / 0.25;
                return [0.1, 1.0, 1.0 - u];
            }} else if (t < 0.75) {{
                let u = (t - 0.5) / 0.25;
                return [u, 1.0, 0.0];
            }} else {{
                let u = (t - 0.75) / 0.25;
                return [1.0, 1.0 - 0.8 * u, 0.0];
            }}
        }}

        const pointCount = POINT_DEVS.length;
        const pointGeometry = new THREE.BufferGeometry();
        const positions = new Float32Array(POINT_COORDS);
        const colors = new Float32Array(pointCount * 3);

        for (let i = 0; i < pointCount; i++) {{
            const dev = POINT_DEVS[i];
            const normDev = (dev + TOLERANCE_BAND) / (2.0 * TOLERANCE_BAND);
            const [r, g, b] = turboColormap(normDev);
            colors[i * 3 + 0] = r;
            colors[i * 3 + 1] = g;
            colors[i * 3 + 2] = b;
        }}

        pointGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        pointGeometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        pointGeometry.computeBoundingSphere();

        const pointMaterial = new THREE.PointsMaterial({{
            size: 2.5,
            vertexColors: true,
            sizeAttenuation: true,
            clippingPlanes: [clipPlaneZ]
        }});

        const pointCloud = new THREE.Points(pointGeometry, pointMaterial);
        scene.add(pointCloud);

        let cadMesh = null;
        if (MESH_VERTS.length > 0 && MESH_TRIS.length > 0) {{
            const meshGeo = new THREE.BufferGeometry();
            meshGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(MESH_VERTS), 3));
            meshGeo.setIndex(MESH_TRIS);
            meshGeo.computeVertexNormals();

            const meshMat = new THREE.MeshStandardMaterial({{
                color: 0x94a3b8,
                roughness: 0.4,
                metalness: 0.2,
                transparent: true,
                opacity: 0.4,
                side: THREE.DoubleSide,
                clippingPlanes: [clipPlaneZ]
            }});

            cadMesh = new THREE.Mesh(meshGeo, meshMat);
            scene.add(cadMesh);
        }}

        if (pointGeometry.boundingSphere) {{
            controls.target.copy(pointGeometry.boundingSphere.center);
            camera.position.set(
                pointGeometry.boundingSphere.center.x + pointGeometry.boundingSphere.radius * 1.8,
                pointGeometry.boundingSphere.center.y + pointGeometry.boundingSphere.radius * 1.4,
                pointGeometry.boundingSphere.center.z + pointGeometry.boundingSphere.radius * 1.8
            );
            controls.update();
        }}

        window.addEventListener('resize', () => {{
            camera.aspect = container.clientWidth / container.clientHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(container.clientWidth, container.clientHeight);
        }});

        const raycaster = new THREE.Raycaster();
        raycaster.params.Points.threshold = 1.5;
        const mouse = new THREE.Vector2();
        const probeCard = document.getElementById('probe-card');

        container.addEventListener('mousemove', (e) => {{
            const rect = container.getBoundingClientRect();
            mouse.x = ((e.clientX - rect.left) / container.clientWidth) * 2 - 1;
            mouse.y = -((e.clientY - rect.top) / container.clientHeight) * 2 + 1;

            raycaster.setFromCamera(mouse, camera);
            const intersects = raycaster.intersectObject(pointCloud);

            if (intersects.length > 0) {{
                const idx = intersects[0].index;
                document.getElementById('probe-x').textContent = POINT_COORDS[idx * 3 + 0].toFixed(3);
                document.getElementById('probe-y').textContent = POINT_COORDS[idx * 3 + 1].toFixed(3);
                document.getElementById('probe-z').textContent = POINT_COORDS[idx * 3 + 2].toFixed(3);
                const dev = POINT_DEVS[idx];
                const probeDev = document.getElementById('probe-dev');
                probeDev.textContent = (dev >= 0 ? "+" : "") + dev.toFixed(4);
                probeDev.style.color = Math.abs(dev) > TOLERANCE_BAND ? "#ef4444" : "#10b981";

                probeCard.style.display = 'block';
                probeCard.style.left = (e.clientX - rect.left + 15) + 'px';
                probeCard.style.top = (e.clientY - rect.top + 15) + 'px';
            }} else {{
                probeCard.style.display = 'none';
            }}
        }});

        document.getElementById('toggle-cad').addEventListener('change', (e) => {{
            if (cadMesh) cadMesh.visible = e.target.checked;
        }});
        document.getElementById('toggle-wireframe').addEventListener('change', (e) => {{
            if (cadMesh) cadMesh.material.wireframe = e.target.checked;
        }});
        document.getElementById('toggle-points').addEventListener('change', (e) => {{
            pointCloud.visible = e.target.checked;
        }});
        document.getElementById('slider-point-size').addEventListener('input', (e) => {{
            pointMaterial.size = parseFloat(e.target.value);
            document.getElementById('point-size-val').textContent = e.target.value + " px";
        }});
        document.getElementById('slider-clip-z').addEventListener('input', (e) => {{
            const val = parseFloat(e.target.value);
            if (val >= 99) {{
                clipPlaneZ.constant = 1000;
                document.getElementById('clip-z-val').textContent = "Off";
            }} else {{
                clipPlaneZ.constant = val;
                document.getElementById('clip-z-val').textContent = val.toFixed(1) + " mm";
            }}
        }});

        function drawHistogram() {{
            const histCanvas = document.getElementById('histogram-canvas');
            const ctx = histCanvas.getContext('2d');
            const w = histCanvas.width = histCanvas.clientWidth * window.devicePixelRatio;
            const h = histCanvas.height = histCanvas.clientHeight * window.devicePixelRatio;

            const numBins = 40;
            const bins = new Array(numBins).fill(0);
            const range = TOLERANCE_BAND * 2.0;

            for (let i = 0; i < POINT_DEVS.length; i++) {{
                const d = POINT_DEVS[i];
                let bin = Math.floor(((d + range) / (range * 2)) * numBins);
                if (bin < 0) bin = 0;
                if (bin >= numBins) bin = numBins - 1;
                bins[bin]++;
            }}

            const maxCount = Math.max(...bins, 1);
            const barWidth = w / numBins;

            ctx.clearRect(0, 0, w, h);
            for (let b = 0; b < numBins; b++) {{
                const barHeight = (bins[b] / maxCount) * (h - 10);
                const x = b * barWidth;
                const y = h - barHeight;

                const centerDev = -range + (b + 0.5) * (range * 2 / numBins);
                if (Math.abs(centerDev) <= TOLERANCE_BAND) {{
                    ctx.fillStyle = '#10b981';
                }} else {{
                    ctx.fillStyle = '#ef4444';
                }}

                ctx.fillRect(x, y, barWidth - 1, barHeight);
            }}

            const zeroX = (range / (range * 2)) * w;
            ctx.strokeStyle = '#ffffff';
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(zeroX, 0);
            ctx.lineTo(zeroX, h);
            ctx.stroke();
        }}
        drawHistogram();

        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }}
        animate();
    </script>
</body>
</html>
"""
    return html


def save_html_viewer(
    result: InspectionResult,
    output_path: str | Path,
    title: Optional[str] = None,
    tolerance_band: float = 0.5,
) -> Path:
    """Save an interactive HTML WebGL inspection dashboard to disk.

    Args:
        result: The InspectionResult object.
        output_path: Path to target .html file.
        title: Optional custom dashboard title.
        tolerance_band: Color normalization range in mm.

    Returns:
        The written Path.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    html_content = generate_html_viewer(result, title=title, tolerance_band=tolerance_band)
    path.write_text(html_content, encoding="utf-8")
    return path


# Convenient alias for top-level export
generate_inspection_dashboard = generate_html_viewer


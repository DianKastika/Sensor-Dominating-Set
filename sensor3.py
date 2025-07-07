import os
import streamlit as st
import geopandas as gpd
import folium
from shapely.geometry import Point
import numpy as np
import pandas as pd
from streamlit_folium import st_folium
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from tempfile import TemporaryDirectory

st.set_page_config(layout="wide")
st.title("Penempatan Sensor dengan Dominating Set")
st.sidebar.header("Pilih Lokasi Penempatan Sensor")

sensor_mode = st.sidebar.radio("Lokasi Sensor", ["Sensor di Laut", "Sensor di Darat"])
grid_spacing = st.sidebar.number_input("Ukuran Grid (derajat)", 0.01, 1.0, 0.1)

# === Parameter Radius Sensor ===
if sensor_mode == "Sensor di Laut":
    radius_val = st.sidebar.number_input("Jangkauan Sensor (nautical miles)", 1, 300, 50)
    min_distance_val = st.sidebar.number_input("Jarak Minimum antar Sensor (nautical miles)", 1, 300, 10)
    radius_deg = radius_val / 60
    min_distance_deg = min_distance_val / 60
    satuan_jarak = "nautical miles"
    radius_km = round(radius_val * 1.852, 1)
    min_distance_km = round(min_distance_val * 1.852, 1)
else:
    radius_val = st.sidebar.number_input("Jangkauan Sensor (km)", 1, 500, 50)
    min_distance_val = st.sidebar.number_input("Jarak Minimum antar Sensor (km)", 1, 500, 10)
    radius_deg = radius_val / 111
    min_distance_deg = min_distance_val / 111
    satuan_jarak = "km"
    radius_km = radius_val
    min_distance_km = min_distance_val

# === Fungsi Validasi dan Loader ===
def load_vector_file(uploaded_file):
    if uploaded_file is None:
        return None
    filename = uploaded_file.name
    try:
        if filename.endswith(".geojson") or filename.endswith(".json"):
            gdf = gpd.read_file(uploaded_file)
        elif filename.endswith(".zip"):
            with TemporaryDirectory() as tmpdir:
                temp_path = os.path.join(tmpdir, filename)
                with open(temp_path, "wb") as f:
                    f.write(uploaded_file.read())
                gdf = gpd.read_file(f"zip://{temp_path}")
        else:
            st.error("⚠️ Format file tidak dikenali. Gunakan .geojson, .json, atau shapefile .zip")
            return None
        if gdf.empty:
            st.error("⚠️ File tidak mengandung geometri.")
            return None
        return gdf
    except Exception as e:
        st.error(f"⚠️ Gagal memuat file: {e}")
        return None

def load_csv_sensors(file):
    sensors = []
    if file:
        try:
            df = pd.read_csv(file)
            if not {'longitude', 'latitude'}.issubset(df.columns):
                st.error("⚠️ CSV harus memiliki kolom 'latitude' dan 'longitude'")
                return []
            for _, row in df.iterrows():
                sensors.append(Point(row['longitude'], row['latitude']))
        except Exception as e:
            st.error(f"⚠️ Gagal membaca CSV: {e}")
    return sensors

def generate_grid_points(minx, miny, maxx, maxy, spacing):
    xs = np.arange(minx, maxx, spacing)
    ys = np.arange(miny, maxy, spacing)
    return gpd.GeoDataFrame(geometry=[Point(x, y) for x in xs for y in ys], crs="EPSG:4326")

# === Inisialisasi
intersection_points = []
sensors = []

# === SENSOR DARAT ===
if sensor_mode == "Sensor di Darat":
    batas_file = st.sidebar.file_uploader("Unggah Batas Wilayah (GeoJSON / Shapefile dalam ZIP)", type=["geojson", "json", "zip"])
    darat_file = st.sidebar.file_uploader("Unggah Wilayah Darat (GeoJSON / Shapefile dalam ZIP)", type=["geojson", "json", "zip"])
    sensor_file = st.sidebar.file_uploader("Sensor Eksisting (CSV, opsional)", type="csv")

    if batas_file and darat_file:
        batas = load_vector_file(batas_file)
        darat = load_vector_file(darat_file)

        if batas is not None and darat is not None:
            batas_union = unary_union(batas.geometry)
            darat_union = unary_union(darat.geometry)
            minx, miny, maxx, maxy = batas.total_bounds
            grid = generate_grid_points(minx, miny, maxx, maxy, grid_spacing)
            grid = grid[grid.geometry.within(darat_union)]
            grid = grid[grid.geometry.within(batas_union)]

            if grid.empty:
                st.warning("⚠️ Tidak ada titik grid yang valid di wilayah darat.")
            else:
                intersection_points = list(grid.geometry)
                sensors = load_csv_sensors(sensor_file)

# === SENSOR LAUT ===
elif sensor_mode == "Sensor di Laut":
    batas_file = st.sidebar.file_uploader("Unggah Batas Wilayah (GeoJSON / ZIP)", type=["geojson", "json", "zip"])
    garis_file = st.sidebar.file_uploader("Unggah Garis Pantai (GeoJSON / ZIP)", type=["geojson", "json", "zip"])
    sensor_file = st.sidebar.file_uploader("Sensor Eksisting (CSV, opsional)", type="csv")

    if batas_file and garis_file:
        batas = load_vector_file(batas_file)
        garis = load_vector_file(garis_file)

        if batas is not None and garis is not None:
            pantai_union = unary_union(garis.geometry)
            batas_union = unary_union(batas.geometry)
            grid = generate_grid_points(*garis.total_bounds, grid_spacing)
            grid = grid[grid.geometry.intersects(pantai_union)]
            grid = grid[grid.geometry.within(batas_union)]

            if grid.empty:
                st.warning("⚠️ Tidak ada titik grid yang memotong garis pantai.")
            else:
                intersection_points = list(grid.geometry)
                sensors = load_csv_sensors(sensor_file)

# === PROSES DOMINASI ===
if intersection_points:
    coords = np.array([[p.x, p.y] for p in intersection_points])
    tree = cKDTree(coords)
    coverage = np.zeros(len(coords), dtype=int)

    for pt in sensors:
        idxs = tree.query_ball_point([pt.x, pt.y], r=radius_deg)
        for i in idxs:
            coverage[i] += 1

    remaining_idx = np.where(coverage == 0)[0]
    remaining_coords = coords[remaining_idx]
    selected_sensors = sensors.copy()

    max_iter = 1000
    for _ in range(max_iter):
        if len(remaining_idx) == 0:
            break
        tree_remain = cKDTree(remaining_coords)
        counts = tree_remain.query_ball_point(remaining_coords, r=radius_deg, return_length=True)
        sorted_idx = np.argsort(-np.array(counts))
        for idx in sorted_idx:
            candidate = Point(remaining_coords[idx])
            if all(candidate.distance(p) >= min_distance_deg for p in selected_sensors):
                idxs = tree.query_ball_point([candidate.x, candidate.y], r=radius_deg)
                new_covered = [i for i in idxs if coverage[i] == 0]
                if not new_covered:
                    continue
                selected_sensors.append(candidate)
                for i in new_covered:
                    coverage[i] += 1
                remaining_idx = np.where(coverage == 0)[0]
                remaining_coords = coords[remaining_idx]
                break

    # === VISUALISASI ===
    m = folium.Map(location=[-2.5, 118], zoom_start=5)
    for pt in sensors:
        folium.Marker(location=[pt.y, pt.x], icon=folium.Icon(color='blue')).add_to(m)
        folium.Circle(location=[pt.y, pt.x], radius=radius_deg * 111000, color='blue', fill=True, fill_opacity=0.1).add_to(m)

    new_sensors = [pt for pt in selected_sensors if pt not in sensors]
    for pt in new_sensors:
        folium.Marker(location=[pt.y, pt.x], icon=folium.Icon(color='orange')).add_to(m)
        folium.Circle(location=[pt.y, pt.x], radius=radius_deg * 111000, color='orange', fill=True, fill_opacity=0.1).add_to(m)

    st.subheader("🗺️ Peta Sensor Keamanan")
    st_folium(m, width=1200, height=600)

    st.subheader("📊 Statistik Penempatan Sensor")
    st.markdown(f"""
    - Total titik potensial: **{len(coords)}**
    - Sensor eksisting: **{len(sensors)}**
    - Sensor baru terpilih: **{len(new_sensors)}**
    - Total sensor setelah dominasi: **{len(selected_sensors)}**
    - Radius jangkauan: **{radius_val} {satuan_jarak}** (~{radius_km} km)
    - Jarak minimum antar sensor: **{min_distance_val} {satuan_jarak}** (~{min_distance_km} km)
    """)

    df_result = pd.DataFrame({
        'latitude': [pt.y for pt in selected_sensors],
        'longitude': [pt.x for pt in selected_sensors],
        'radius': [radius_val] * len(selected_sensors),
        'satuan_radius': [satuan_jarak] * len(selected_sensors)
    })
    st.download_button("📥 Unduh CSV Titik Sensor", data=df_result.to_csv(index=False).encode(), file_name="sensor_output.csv", mime="text/csv")
    html = m.get_root().render()
    st.download_button("🌐 Unduh Peta Interaktif (HTML)", data=html.encode(), file_name="peta_sensor_interaktif.html", mime="text/html")
else:
    st.warning("⚠️ Silakan unggah file dan pastikan titik potensial berhasil dibuat.")


import matplotlib.pyplot as plt
import pandas as pd
from bs4 import BeautifulSoup
import requests
import re
from urllib.parse import urljoin
import numpy as np
from io import StringIO
import time
import sqlite3
from tqdm import tqdm


patron = r"(.+?)\s+(\d{2}:\d{2}:\d+\.\d+)\s+([+-]\d{2}:\d{2}:\d+\.\d+)"

# Paso 1: Abrir el archivo HTML y leer su contenido
def leer_ccd(numero_ccd):
    datos = []
    url = f"http://rubin.astro.puc.cl:3672/dosc/Shallow24A/3672/ccd{numero_ccd}.html"
    respuesta = requests.get(url)

    if respuesta.status_code != 200: #la 200 es la buena

        return None
    
    soup = BeautifulSoup(respuesta.text, "html.parser")

    # Codigo para extraer datos

    columnas = soup.find_all("div", class_ = "column")
    
    candidatos_tot = len(columnas)

    for candidato in tqdm(columnas, desc=f"Leyendo CCD {numero_ccd}"):
        # Estos son todos los candidatos del mismo ccd.
        texto = candidato.find("div", class_="text-overlay").get_text(" ", strip=True) 
        m = re.search(patron, texto)
        links = candidato.find_all("a")

        
        if m:
            nombre = m.group(1)
            categoria = None
            ra = m.group(2)
            dec = m.group(3) 
            links = candidato.find_all("a")
            dosc_id = f"{numero_ccd:02d}_{ra.replace(':','')}{dec.replace(':','')}"
            url_candidato = None

            # Para encontrar el url del candidato
            for link in links:
                href = link.get("href")
                if "href" and "index.html" in href:
                    url_candidato = href
                    break
            ned_source =  None #leer_NED_source(url_candidato)
            
            
            nearest_type, nearest_offset, nearest_filter = obtener_tipo_nearest(candidato)
            
            m2 = re.search(r"(.+?)\s*\((.+)\)", nombre)
            if m2:
                nombre = m2.group(1).strip()
                categoria = m2.group(2).strip()
            

            datos.append({
                "photpipe_id": nombre,
                "ccd": numero_ccd,
                "categoria": categoria,
                "ra": ra,
                "dec": dec,
                "url": url_candidato,
                "dosc_id": dosc_id,
                "nearest_type": nearest_type,
                "ned_source": ned_source
            })



    return pd.DataFrame(datos)

def leer_forcedlc(url_forcedlc):
    respuesta = requests.get(url_forcedlc)
    

    df = pd.read_csv(
    url_forcedlc,
    sep=r"\s+",
    comment="#",
    names=[
        "MJD",
        "dateobs",
        "photcode",
        "filt",
        "flux_c",
        "dflux_c",
        "type",
        "chisqr",
        "ZPTMAG_c",
        "m",
        "dm",
        "ra",
        "dec",
        "cmpfile",
        "tmpl"
    ],
    engine="python"
                    )

    return df


def obtener_forcedlc(url_candidato, photpipe_id):
    """"
    Esta funcion lo que hace es encontrar la url al forcedlc. Solo eso. 
    No extrae ningún tipo de informacion o data más que eso.
    """
    respuesta = requests.get(url_candidato)
    if respuesta.status_code != 200:
        print("No pude entrar al candidato")
        return None
    
    # Extraer numero del ID
    candidato_id = photpipe_id.split("_")[-1] 
    
    soup = BeautifulSoup(respuesta.text, "html.parser")

    enlaces = soup.find_all("a")
    
    for enlace in enlaces:
        href = enlace.get("href")
        if href and f"cand{candidato_id}.forced.difflc.txt" in href:
            return href
    
    return None

def graficar_forcedlc(candidato):
    
    colores = {
        "g": "blue",
        "r": "green",
        "i": "red",
        "z": "purple"
    }

    marcadores = {
        "g": "D",
        "r": "o",
        "i": "s",
        "z": "^"
    }


    # ----------------------------
    # Obtener forced light curve
    # ----------------------------

    url_candidato = candidato["url"]
    photpipe_id = candidato["photpipe ID"]

    href_forcedlc = obtener_forcedlc(url_candidato, photpipe_id)

    if href_forcedlc is None:
        print("No se encontró forcedlc")
        return

    url_forcedlc = urljoin(url_candidato, href_forcedlc)
    df = leer_forcedlc(url_forcedlc)


    # ----------------------------
    # Preparar datos
    # ----------------------------

    df["m"] = pd.to_numeric(df["m"], errors="coerce")
    df["dm"] = pd.to_numeric(df["dm"], errors="coerce")

    df["m_calculada"] = (
        df["ZPTMAG_c"] - 2.5*np.log10(df["flux_c"])
    )

    df["mag_final"] = df["m"]

    df.loc[
        df["mag_final"].isna(),
        "mag_final"
    ] = df.loc[
        df["mag_final"].isna(),
        "m_calculada"
    ]


    # solo magnitudes con error
    df_mag = df.dropna(subset=["m", "dm"])
    peak_mjd, peak_mag, peak_filtro = encontrar_peak_magnitud(df)

    # ----------------------------
    # Crear figura
    # ----------------------------

    fig, axes = plt.subplots(
        2, 1,
        figsize=(10,8),
        sharex=True
    )


    # ----------------------------
    # Magnitud
    # ----------------------------

    for banda, color in colores.items():

        datos = df_mag[df_mag["filt"] == banda]

        axes[0].errorbar(
            datos["MJD"],
            datos["mag_final"],
            yerr=datos["dm"],
            fmt=marcadores[banda],
            markersize=6,
            markerfacecolor="none",
            markeredgecolor=color,
            color=color,
            label=banda
        )


    axes[0].invert_yaxis()
    axes[0].set_ylabel("Magnitude")
    axes[0].set_ylim(25.5,22)
    axes[0].set_title(
        f"Forced light curve - {photpipe_id}"
    )
    axes[0].legend()
    axes[0].grid(True)
    axes[0].scatter(
        peak_mjd,
        peak_mag,
        marker="*",
        s=150,
        color="black",
        label="Peak"
                    )

    axes[0].legend()


    # ----------------------------
    # Flux
    # ----------------------------

    for banda, color in colores.items():

        datos = df[df["filt"] == banda]

        axes[1].plot(
            datos["MJD"],
            datos["flux_c"],
            linestyle="None",
            marker=marcadores[banda],
            markersize=7,
            markerfacecolor="none",
            markeredgecolor=color,
            color=color,
            label=banda
        )


    axes[1].set_xlabel("MJD")
    axes[1].set_ylabel("Flux")
    axes[1].set_title("Forced light curve - Flux")
    axes[1].grid(True)

    axes[1].axhline(
        y=0,
        color="gray",
        linewidth=2,
        linestyle=":"
    )


    plt.tight_layout()
    plt.show()

def obtener_dosc_id(candidato):
    
    ccd = int(candidato["CCD"])
    
    ra = candidato["RA"].replace(":", "")
    dec = candidato["DEC"].replace(":", "")
    
    dosc_id = f"{ccd:02d}_{ra}{dec}"
    
    return dosc_id

def encontrar_peak_magnitud(df):
    
    # Solo datos con magnitud válida
    datos = df.dropna(subset=["mag_final"])

    # encontrar índice de menor magnitud
    idx_peak = datos["mag_final"].idxmin()

    peak_mjd = datos.loc[idx_peak, "MJD"]
    peak_mag = datos.loc[idx_peak, "mag_final"]
    peak_filtro = datos.loc[idx_peak, "filt"]

    return peak_mjd, peak_mag, peak_filtro

def obtener_tipo_nearest(candidato):

    texto = candidato.get_text(" ", strip=True)

    patron = r"([griz])=([\d.]+)\"\s+\((Star|Gal)\)"

    matches = re.findall(patron, texto)

    if matches:

        offsets = []

        for filtro, offset, tipo in matches:
            offsets.append({
                "filtro": filtro,
                "offset": float(offset),
                "tipo": tipo
            })

        # buscar el menor offset
        minimo = min(offsets, key=lambda x: x["offset"])

        return minimo["tipo"], minimo["offset"], minimo["filtro"]

    return None, None, None

def obtener_NED_nearest_source(url_candidato):
    
    respuesta = requests.get(url_candidato)

    soup = BeautifulSoup(respuesta.text, "html.parser")

    links = soup.find_all("a")

    url_ned = None

    for link in links:
        href = link.get("href")
        if href and "ned.ipac.caltech.edu" in href:
            url_ned = href
            break

    return url_ned

def leer_NED_source(url_ned):

    respuesta = requests.get(url_ned)

    soup = BeautifulSoup(respuesta.text, "html.parser")

    tablas = pd.read_html(StringIO(respuesta.text))


def guardar_candidatos_sqlite(df):

    conn = sqlite3.connect("supernovae_v2.db")
    cursor = conn.cursor()

    agregados = 0
    duplicados = 0

    for _, fila in tqdm(
        df.iterrows(),
        total=len(df),
        desc="Guardando candidatos"
    ):

        cursor.execute("""
        INSERT OR IGNORE INTO candidates
        (
        photpipe_id,
        ccd,
        categoria,
        ra,
        dec,
        url,
        dosc_id,
        nearest_type,
        ned_source
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
        fila["photpipe_id"],
        fila["ccd"],
        fila["categoria"],
        fila["ra"],
        fila["dec"],
        fila["url"],
        fila["dosc_id"],
        fila.get("nearest_type"),
        fila.get("ned_source")
        ))


        if cursor.rowcount == 1:
            agregados += 1
        else:
            duplicados += 1


    conn.commit()
    conn.close()


    return {
        "agregados": agregados,
        "duplicados": duplicados
    }
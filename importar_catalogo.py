from funciones_sn import (leer_ccd, 
                          guardar_candidatos_sqlite,
                          )  
import time 


def cargar_ccd(numero_ccd):
    
    print(f"\nCargando CCD {numero_ccd}...")

    df = leer_ccd(numero_ccd)
    

    if df is None or len(df) == 0:
        print("No se encontraron candidatos")
        return


    print(f"Candidatos encontrados: {len(df)}")


    resultado = guardar_candidatos_sqlite(df)

    print("\nCarga terminada")
    print(f"Agregados: {resultado['agregados']}")
    print(f"Ya existían: {resultado['duplicados']}")

if __name__ == "__main__":


    for i in range(1, 61):
        cargar_ccd(i)
import numpy as np
from read_gff_header_1_8b import read_gff_header_1_8b


def load_gff_1_8b(FileName: str, amplitude_only: bool = True):
    """
    Загружает GFF-изображение.

    Параметры
    ----------
    FileName : str
        Путь к файлу.
    amplitude_only : bool, optional
        Если True, возвращает амплитуду (sqrt(magnitude)), иначе комплексное число.

    Возвращает
    ----------
    Image : np.ndarray | None
        Массив пикселей или None при ошибке.
    Header : dict | None
        Заголовок файла или None при ошибке.
    """
    # ---- Открытие файла и чтение заголовка ----
    Header, fid = read_gff_header_1_8b(FileName)

    # ---- Проверка ошибок ----
    if isinstance(fid, int) and fid < 0:
        return None, None

    # ---- Порядок индексов ----
    if Header["RowMajor"]:
        FirstIndexLength = Header["RgCnt"]
        SecondIndexLength = Header["AzCnt"]
        Flip = False
    else:
        FirstIndexLength = Header["AzCnt"]
        SecondIndexLength = Header["RgCnt"]
        Flip = True

    Image = None

    # ---- Разбор по типу изображения ----
    if Header["ImageType"] == 0:
        print("Reading sqrt(magnitude) pixels.")
        n = FirstIndexLength * SecondIndexLength
        raw = fid.read(n)
        Image = np.frombuffer(raw, dtype=np.uint8)

    elif Header["ImageType"] == 1:
        bpp = int(Header["BytesPerPixel"])
        if bpp == 4:
            # phase = 2 байта (uint16), mag = 2 байта (uint16)
            n = FirstIndexLength * 2 * SecondIndexLength
            raw = fid.read(n * 2)   # 2 байта на ushort
            arr = np.frombuffer(raw, dtype="<u2" if Header["Endian"] == 0 else ">u2")
            if amplitude_only:
                print("Reading sqrt(magnitude) pixels.")
                Image = np.sqrt(arr[1::2].astype(np.float32))
            else:
                print("Reading complex (magnitude and phase) pixels.")
                phase = arr[0::2].astype(np.float64)
                mag = arr[1::2].astype(np.float64)
                Image = (mag * np.exp(1j * phase * 2 * np.pi / (2 ** 16))).astype(np.complex64)
        elif bpp == 8:
            # phase = 4 байта (uint32), mag = 4 байта (uint32)
            n = FirstIndexLength * 2 * SecondIndexLength
            raw = fid.read(n * 4)
            arr = np.frombuffer(raw, dtype="<u4" if Header["Endian"] == 0 else ">u4")
            if amplitude_only:
                print("Reading sqrt(magnitude) pixels.")
                Image = np.sqrt(arr[1::2].astype(np.float32))
            else:
                print("Reading complex (magnitude and phase) pixels.")
                phase = arr[0::2].astype(np.float64)
                mag = arr[1::2].astype(np.float64)
                Image = (mag * np.exp(1j * phase * 2 * np.pi / (2 ** 32))).astype(np.complex64)
        else:
            raise ValueError(
                f"ImageType=1: неподдерживаемый BytesPerPixel={bpp} "
                f"(ожидается 4 или 8)"
            )

    elif Header["ImageType"] == 2:
        n = 2 * FirstIndexLength * SecondIndexLength
        raw = fid.read(n * 4)   # float32
        arr = np.frombuffer(raw, dtype="<f4" if Header["Endian"] == 0 else ">f4")
        if amplitude_only:
            print("Reading sqrt(real) pixels.")
            Image = np.sqrt(arr[0::2].astype(np.float32))
        else:
            print("Reading complex (real and imaginary) pixels.")
            real = arr[0::2].astype(np.float32)
            imag = arr[1::2].astype(np.float32)
            Image = (real + 1j * imag).astype(np.complex64)
    else:
        print(f"Unknown pixel type: {Header['ImageType']}.", file=sys.stderr)

    if Image is not None:
        Image = Image.reshape(SecondIndexLength, FirstIndexLength)
        if Flip:
            Image = Image.T

    fid.close()

    return Image, Header

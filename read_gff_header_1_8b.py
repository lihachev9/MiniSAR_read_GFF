import os
import struct
import sys


def _read(fid, fmt: str, size: int):
    """Прочитать size байт и распаковать по формату fmt."""
    data = fid.read(size)
    if len(data) != size:
        raise EOFError(f"Недостаточно данных: ожидалось {size}, получено {len(data)}")
    return struct.unpack(fmt, data)


def _read_char(fid, n: int) -> str:
    """Прочитать n байт как строку (обрезая нули)."""
    data = fid.read(n)
    return data.rstrip(b"\x00").decode("latin-1")


def read_gff_header_1_8b(FileName: str):
    """
    Читает заголовок GFF-файла.

    Параметры
    ----------
    FileName : str
        Имя файла.

    Возвращает
    ----------
    Header : dict | None
        Словарь с полями заголовка. None при ошибке.
    fid_out : file | int
        Открытый файл (позиционирован на начало пиксельных данных)
        либо код ошибки:
          -1 — файл не удалось открыть,
          -2 — файл не существует,
          -3 — файл не является GFF,
    """

    if not os.path.isfile(FileName):
        print(f"GSAT File {FileName} does not exist.", file=sys.stderr)
        return None, -2
    try:
        fid = open(FileName, "rb")
    except OSError:
        print(f"GSAT File {FileName} cannot be opened.", file=sys.stderr)
        return None, -1

    # ---- Проверка "GSATIMG" ----
    fid.seek(0)
    struct_id = _read_char(fid, 7)
    if struct_id != "GSATIMG":
        print("Selected file is not in GFF.", file=sys.stderr)
        fid.close()
        return None, -3

    # ---- Определение endianness ----
    fid.seek(54)
    (endian_flag,) = _read(fid, "<H", 2)  # 0 = little-endian
    fid.close()

    if endian_flag == 0:
        endian = "<"
        print(f"Opened {FileName} as ieee-le")
    else:
        endian = ">"
        print(f"Opened {FileName} as ieee-be")

    fid = open(FileName, "rb")
    Header = {}

    # ---- Версия ----
    fid.seek(8)
    Header["Version_Minor"] = _read(fid, endian + "H", 2)[0]
    Header["Version_Major"] = _read(fid, endian + "H", 2)[0]
    Header["Length"] = _read(fid, endian + "I", 4)[0]

    Header["Creator_Length"] = _read(fid, endian + "H", 2)[0]
    Header["Creator"] = _read_char(fid, Header["Creator_Length"])
    fid.read(24 - Header["Creator_Length"])

    if fid.tell() != 42:
        print("file marker not at 42 to start date & time.", file=sys.stderr)

    Header["DateTime"] = list(_read(fid, endian + "6H", 12))  # yr, mo, da, hr, min, sec
    Header["Endian"] = _read(fid, endian + "H", 2)[0]

    vmajor = Header["Version_Major"]
    vminor = Header["Version_Minor"]

    if (vmajor == 1 and vminor > 7) or vmajor > 1:
        Header["BytesPerPixel"] = _read(fid, endian + "f", 4)[0]
    else:
        Header["BytesPerPixel"] = _read(fid, endian + "I", 4)[0]

    Header["FrameCnt"] = _read(fid, endian + "I", 4)[0]
    Header["ImageType"] = _read(fid, endian + "I", 4)[0]
    Header["RowMajor"] = _read(fid, endian + "I", 4)[0]
    Header["RgCnt"] = _read(fid, endian + "I", 4)[0]
    Header["AzCnt"] = _read(fid, endian + "I", 4)[0]
    Header["ScaleExponent"] = _read(fid, endian + "i", 4)[0]
    Header["ScaleMantissa"] = _read(fid, endian + "i", 4)[0]
    Header["OffsetExponent"] = _read(fid, endian + "i", 4)[0]
    Header["OffsetMantissa"] = _read(fid, endian + "i", 4)[0]
    fid.read(32)  # Res2 filler

    if fid.tell() != 128:
        print("file marker not at 128 to start comment.", file=sys.stderr)

    Header["Comment_Length"] = _read(fid, endian + "H", 2)[0]
    Header["Comment"] = _read_char(fid, Header["Comment_Length"])
    fid.read(166 - Header["Comment_Length"])

    Header["ImagePlane"] = _read(fid, endian + "I", 4)[0]
    Header["RgPixelSz"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["AzPixelSz"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["AzOverlap"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-16)
    Header["SRPLat"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["SRPLong"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["SRPAlt"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-16)
    Header["RFOA"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["XtoSRP"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-16)
    fid.read(32)  # Res3 filler

    if fid.tell() != 364:
        print("file marker not at 364 to start phase history ID.", file=sys.stderr)

    Header["PhaseName_Length"] = _read(fid, endian + "H", 2)[0]
    Header["PhaseName"] = _read_char(fid, Header["PhaseName_Length"]).strip()
    fid.read(128 - Header["PhaseName_Length"])
    Header["ImageName_Length"] = _read(fid, endian + "H", 2)[0]
    Header["ImageName"] = _read_char(fid, Header["ImageName_Length"]).strip()
    fid.read(128 - Header["ImageName_Length"])

    Header["LookCnt"] = _read(fid, endian + "I", 4)[0]
    Header["ParamRefAp"] = _read(fid, endian + "I", 4)[0]
    Header["ParamRefPos"] = _read(fid, endian + "I", 4)[0]
    Header["GrzAngle"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-23)
    Header["Squint"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["GTA"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["RgBeamCtr"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-8)
    Header["FlightTime"] = _read(fid, endian + "I", 4)[0] * 10.0 ** (-3)

    if (vmajor == 1 and vminor > 5) or vmajor > 1:
        Header["RgChirpRate"] = _read(fid, endian + "f", 4)[0]
    else:
        Header["RgChirpRate"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (+12)

    Header["XtoStart"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-16)
    Header["MoCompMode"] = _read(fid, endian + "I", 4)[0]
    Header["V_X"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["APCLat"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["APCLong"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["APCAlt"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["Calparm"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-24)
    Header["LogicalBlkAdrr"] = _read(fid, endian + "I", 4)[0]

    if fid.tell() != 692:
        print("file marker not at 692 to start Azimuth Resolution.", file=sys.stderr)

    Header["AzResolution"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["RgResolution"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
    Header["DesSigmaN"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["DesGrazAngle"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-23)
    Header["DesSquint"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)
    Header["DesRng"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-8)
    Header["SceneTrackAngle"] = _read(fid, endian + "i", 4)[0] * 2.0 ** (-23)

    if fid.tell() != 720:
        print("file marker not at 720 to start User Specified Data.", file=sys.stderr)

    Header["UserParam"] = _read_char(fid, 48)  # User Specified Data

    if fid.tell() != 768:
        print("file marker not at 768 to start CCD Parameters.", file=sys.stderr)

    # CCD Parameters
    for name in (
        "CoarseSNR", "CoarseAzSubSamp", "CoarseRngSubSamp", "MaxAzShift",
        "MaxRngShift", "CoarseDltAz", "CoarseDltRng", "TotProcs",
        "TptBoxCMode", "SNRThresh", "RngSize", "MapBoxSize", "BoxSize",
        "BoxSpc", "TotTPts", "GoodTPts", "RndSeed", "RngShift", "AzShift",
        "SumXRamp", "SumYRamp",
    ):
        Header[name] = _read(fid, endian + "i", 4)[0]

    if fid.tell() != 852:
        print("file marker not at 852 to start General Fields for RTV.", file=sys.stderr)

    Header["Cy9kTapeBlock"] = _read(fid, endian + "I", 4)[0]
    Header["NominalCenterFreq"] = _read(fid, endian + "f", 4)[0]
    Header["ImageFlags"] = _read(fid, endian + "I", 4)[0]
    Header["LineNumber"] = _read(fid, endian + "I", 4)[0]
    Header["PatchNumber"] = _read(fid, endian + "I", 4)[0]
    Header["Lambda0"] = _read(fid, endian + "f", 4)[0]
    Header["SRngPixSpace"] = _read(fid, endian + "f", 4)[0]
    Header["DoppPixSpace"] = _read(fid, endian + "f", 4)[0]
    Header["DoppOffset"] = _read(fid, endian + "f", 4)[0]
    Header["DoppRngScale"] = _read(fid, endian + "f", 4)[0]
    Header["MuxTimeDelay"] = _read(fid, endian + "f", 4)[0]
    Header["APCXECEF"] = _read(fid, endian + "d", 8)[0]
    Header["APCYECEF"] = _read(fid, endian + "d", 8)[0]
    Header["APCZECEF"] = _read(fid, endian + "d", 8)[0]
    Header["VxECEF"] = _read(fid, endian + "f", 4)[0]
    Header["VyECEF"] = _read(fid, endian + "f", 4)[0]
    Header["VzECEF"] = _read(fid, endian + "f", 4)[0]

    if fid.tell() != 932:
        print("file marker not at 932 to start Phase Calibration for IFSAR.", file=sys.stderr)

    Header["PhaseCal"] = _read(fid, endian + "f", 4)[0]
    Header["SRPxECEF"] = _read(fid, endian + "d", 8)[0]
    Header["SRPyECEF"] = _read(fid, endian + "d", 8)[0]
    Header["SRPzECEF"] = _read(fid, endian + "d", 8)[0]
    Header["Res5"] = _read_char(fid, 64)

    # ---- REV 1.8 изменения ----
    if (vmajor == 1 and vminor > 7) or vmajor > 1:
        Header["HeaderLen1"] = _read(fid, endian + "I", 4)[0]
        Header["ImgDate"] = {
            "Year": _read(fid, endian + "H", 2)[0],
            "Month": _read(fid, endian + "H", 2)[0],
            "Day": _read(fid, endian + "H", 2)[0],
            "Hour": _read(fid, endian + "H", 2)[0],
            "Minute": _read(fid, endian + "H", 2)[0],
            "Second": _read(fid, endian + "H", 2)[0],
        }
        Header["CompFileName"] = _read_char(fid, 128)
        Header["RefFileName"] = _read_char(fid, 128)
        Header["IEPlatform"] = _read_char(fid, 24)
        Header["IEProcID"] = _read_char(fid, 12)
        Header["IERadarModel"] = _read_char(fid, 12)
        Header["IERadarID"] = _read(fid, endian + "I", 4)[0]
        Header["IESWID"] = _read_char(fid, 24)
        Header["IFPlatform"] = _read_char(fid, 24)
        Header["IFProcID"] = _read_char(fid, 12)
        Header["IFRadarModel"] = _read_char(fid, 12)
        Header["IFRadarID"] = _read(fid, endian + "I", 4)[0]
        Header["IFSWID"] = _read_char(fid, 24)
        Header["IFAlgo"] = _read_char(fid, 8)
        Header["PHPlatform"] = _read_char(fid, 24)
        Header["PHProcID"] = _read_char(fid, 12)
        Header["PHRadarModel"] = _read_char(fid, 12)
        Header["PHRadarID"] = _read(fid, endian + "I", 4)[0]
        Header["PHSWID"] = _read_char(fid, 24)
        Header["PHDataRcd"] = _read(fid, endian + "I", 4)[0]
        Header["ProcProduct"] = _read(fid, endian + "I", 4)[0]
        Header["MissionText"] = _read_char(fid, 8)
        Header["PHSource"] = _read(fid, endian + "I", 4)[0]
        Header["GPSWeek"] = _read(fid, endian + "I", 4)[0]
        Header["DataCollectReqH"] = _read_char(fid, 14)
        Header["Res6"] = _read_char(fid, 2)
        Header["GridName"] = _read_char(fid, 24)
        Header["PixValLinearity"] = _read(fid, endian + "I", 4)[0]
        Header["ComplexOrReal"] = _read(fid, endian + "I", 4)[0]
        Header["BitsPerMagnitude"] = _read(fid, endian + "H", 2)[0]
        Header["BitsPerPhase"] = _read(fid, endian + "H", 2)[0]
        Header["ComplexOrderType"] = _read(fid, endian + "I", 4)[0]
        Header["PixDataType"] = _read(fid, endian + "I", 4)[0]
        Header["ImageLength"] = _read(fid, endian + "I", 4)[0]
        Header["ImageCmpScheme"] = _read(fid, endian + "I", 4)[0]
        Header["APBO"] = _read(fid, endian + "f", 4)[0]
        Header["AsaPitch"] = _read(fid, endian + "f", 4)[0]
        Header["AsaSquint"] = _read(fid, endian + "f", 4)[0]
        Header["DsaPitch"] = _read(fid, endian + "f", 4)[0]
        Header["IRA"] = _read(fid, endian + "f", 4)[0]
        Header["RxPolarization"] = list(_read(fid, endian + "2f", 8))
        Header["TxPolarization"] = list(_read(fid, endian + "2f", 8))
        Header["VxAvg"] = _read(fid, endian + "f", 4)[0]
        Header["VyAvg"] = _read(fid, endian + "f", 4)[0]
        Header["VzAvg"] = _read(fid, endian + "f", 4)[0]
        Header["APCxAvg"] = _read(fid, endian + "f", 4)[0]
        Header["APCyAvg"] = _read(fid, endian + "f", 4)[0]
        Header["APCzAvg"] = _read(fid, endian + "f", 4)[0]
        Header["AveragingTime"] = _read(fid, endian + "f", 4)[0]
        Header["Dgta"] = _read(fid, endian + "f", 4)[0]
        Header["VelocY"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
        Header["VelocZ"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-16)
        Header["Ba"] = _read(fid, endian + "f", 4)[0]
        Header["Be"] = _read(fid, endian + "f", 4)[0]
        Header["AzGeomCorr"] = _read(fid, endian + "I", 4)[0]
        Header["RngGeomCorr"] = _read(fid, endian + "I", 4)[0]
        Header["AzWinFacBW"] = _read(fid, endian + "f", 4)[0]
        Header["RngWinFacBW"] = _read(fid, endian + "f", 4)[0]
        Header["AzWinID"] = _read_char(fid, 48)
        Header["RngWinID"] = _read_char(fid, 48)
        Header["KeepOutViolPrcnt"] = _read(fid, endian + "f", 4)[0]
        Header["AzCoeff"] = list(_read(fid, endian + "6f", 24))
        Header["PosUncertDown"] = _read(fid, endian + "f", 4)[0]
        Header["PosUncertE"] = _read(fid, endian + "f", 4)[0]
        Header["PosUncertN"] = _read(fid, endian + "f", 4)[0]
        Header["NavAidingType"] = _read(fid, endian + "I", 4)[0]
        Header["TwoDNLPhaseCoeffs"] = list(_read(fid, endian + "10f", 40))
        Header["ClutterSNRThresh"] = _read(fid, endian + "f", 4)[0]
        Header["ElevationCoeff"] = list(_read(fid, endian + "9f", 36))
        Header["MonopulseCoeff"] = list(_read(fid, endian + "12f", 48))
        Header["TwistPntErrPrcnt"] = _read(fid, endian + "f", 4)[0]
        Header["TiltPntErrPrcnt"] = _read(fid, endian + "f", 4)[0]
        Header["AzPntErrPrcnt"] = _read(fid, endian + "f", 4)[0]
        Header["SigmaN"] = _read(fid, endian + "I", 4)[0] * 2.0 ** (-23)
        Header["TakeNum"] = _read(fid, endian + "I", 4)[0]
        Header["IFSARFlags"] = list(_read(fid, endian + "5I", 20))
        Header["MuThreshold"] = _read(fid, endian + "f", 4)[0]
        Header["GffAppType"] = _read(fid, endian + "I", 4)[0]
        Header["Res7"] = _read_char(fid, 8)

    # ---- Позиционирование в конец заголовка ----
    if fid.tell() != Header["Length"]:
        print("File marker not at header length - moving to end.", file=sys.stderr)
        fid.seek(Header["Length"], 0)

    return Header, fid

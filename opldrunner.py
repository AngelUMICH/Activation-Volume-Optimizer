class _OPLDRunner:

    def __init__(self):
        self.opldParameters = {
            "REST": 1,
            "TJOB": 2,
            "TDAM": 2,
            "TPBC": 1,
            "TCRY": 2,
            "NTSTEP": 1,
            "NFREQ": 1,
            "NDEMAX": 5000,
            "NAVMAX": 10000,
            "NDIMER": 60,
            "TLAT": 1,
            "RLAT": 2.8553,
            "RBOX": 30.0,
            "PCUTR": 5.3,
            "DCUTR": 4.21,
            "DECUT": 0.26,
            "NDSMAX": 800,
            "NDRMAX": 3,
            "RDIMER": 0.001,
            "FNRMAX": 0.100,
            "FNRMIN": 0.005,
            "DSTPEF": 0.010,
            "DSTEPM": 0.070,
            "TEMP": 723.000,
            "PKAE": 0.0,
            "DRATE": 0.0001,
            "NFREKL": 300,
            "TIME": 0.000e00,
            "SIGRELON": "true",
            "SIGX": 0.0,
            "SIGY": 0.0,
            "SIGZ": 0.0,
        }

    def __modifyOPLDParameters(self):
        # Modify the parameters for the OPLD code
        self.opldParameters["RLAT"] = self.latticeParameter
        self.opldParameters["RBOX"] = self.crystalCenterPosition[0] * 2

    def createOPLDInputFile(self):

        self.__modifyOPLDParameters()

        with open(f"./par_init.dat", "w") as f:
            for key, value in self.opldParameters.items():
                f.write(f"{value:<14} #{key}\n")

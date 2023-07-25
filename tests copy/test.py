import dpdata



def test_load():
    s = dpdata.LabeledSystem(r"dpdata\tests copy\input.out", fmt= "dftb_plus")
    print(s)
    assert len(s) == 1

if __name__ == "__main__":
    test_load()


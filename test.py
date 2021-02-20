import wgs84conv


def main():
    ecef = wgs84conv.lla2ecef([0.429407141487312, 2.11220668133798, 0]).flatten()
    lla = wgs84conv.ecef2lla(ecef).flatten()

    print(ecef)
    print(lla)


if __name__ == "__main__":
    main()

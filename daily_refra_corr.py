import solar_realtime_pipeline as srp
import argparse
import datetime

yesterday_date_str_in_utc = (datetime.datetime.now(datetime.timezone.utc)
                     - datetime.timedelta(days=1)).strftime('%Y-%m-%d')

def main():
    parser = argparse.ArgumentParser(description='Calculate the daily refraction correction for the LWA')
    parser.add_argument('date', type=str
                        , help='Date to calculate the refraction correction for in the format YYYY-MM-DD'
                        , default=yesterday_date_str_in_utc)

    args = parser.parse_args()

    srp.daily_refra_correction(args.date)

if __name__ == '__main__':
    main()

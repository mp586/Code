DP_months = precipitation_month_avg - precipitation_month_avg_ctl
DP_mo_ROI_NH = DP_months.sel(lat=slice(0.,8.)).sel(lon=slice(0., 10.))
DP_mo_NA = DP_months.sel(lat=slice(15.,25.)).sel(lon=slice(45., 75.))
plt.plot(DP_months.month, DP_mo_ROI_NH, DP_months.month, DP_mo_NA)
plt.plot(DP_months.month, DP_mo_ROI_NH.mean('time'), DP_months.month, DP_mo_NA.mean('time'))
DP_mo_NA
plt.plot(DP_months.month, DP_mo_ROI_NH.mean('lat').mean('lon'), DP_months.month, DP_mo_NA.mean('lat').mean('lon'))
plt.show()
plt.plot(DP_months.month, DP_mo_ROI_NH.mean('lat').mean('lon'), 'k', DP_months.month, DP_mo_NA.mean('lat').mean('lon'), 'r')
plt.show()
%hist

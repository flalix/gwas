// Setup jQuery validator to use Bootstrap-friendly defaults
$.validator.setDefaults({
  errorElement: 'span',
  errorClass: 'help-inline',
  highlight: function (element, errorClass) {
    $(element).parents('div.control-group').addClass('error');
  },
  unhighlight: function (element, errorClass) {
    $(element).parents('div.control-group').removeClass('error');
  }
});

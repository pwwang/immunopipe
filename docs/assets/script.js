(function($) {
    $('.rst-content a[href^="http"]').attr('target', '_blank');
    // Open the FAQ section by default
    // Close them if there are more in the future
    $('.rst-content #faq ~ details').attr('open', true);
})(jQuery);
